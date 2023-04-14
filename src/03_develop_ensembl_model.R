library(GAMBLR)
library(caret)
library(glmnet)
library(tidyverse)
library(vroom)
library(h2o)

###########################################################
##### Define helper functions #############################
###########################################################

# Plot pretty ROC curve
plot_roc <- function(
        h2o_rf_model, # RF model build with h2o
        test_data, # data to use for ROC
        # name of the column containing cluster label
        response_column = "cluster_name",
        negref # reference label
    ){

    # training
    XXX <- h2o.performance(h2o_rf_model)
    metrics <- XXX@metrics[["thresholds_and_metric_scores"]]

    train_data <- data.frame(
        Specificity = metrics$fpr,
        Sensitivity = metrics$tpr,
        data = "Training"
    )

    # testing
    test_data[,response_column] <- ifelse(
        test_data[,response_column]==negref,
        negref,
        paste0("non", negref))

    test_data <- as.h2o(test_data)

    YYY <- h2o.performance(
        h2o_rf_model,
        newdata = test_data
    )

    metrics <- YYY@metrics[["thresholds_and_metric_scores"]]

    test_data <- data.frame(
        Specificity = metrics$fpr,
        Sensitivity = metrics$tpr,
        data = "Testing"
    )

    all_data <- rbind(
        train_data,
        test_data
    )

    # generate plot
    my_roc_plot <- ggplot(
        aes(
            x = Specificity,
            y = Sensitivity,
            color = data
        ),
        data = all_data
      ) +
      geom_line() +
      geom_abline(
          intercept = 0,
          slope = 1,
          size = 0.5,
          linetype = "dashed",
          color = "lightgray"
      ) +
      geom_text(
          x = 0.8,
          y = 0.2,
          label = paste(
            "AUC train:",
            round(XXX@metrics$AUC, 3)
          ),
          color = "blue",
          size = 5,
          family = "Arial"
      ) +
      geom_text(
          x = 0.8,
          y = 0.1,
          label = paste(
            "AUC test:",
            round(YYY@metrics$AUC, 3)
          ),
          color = "blue",
          size = 5,
          family = "Arial"
      ) +
      theme_Morons() +
      ylim(0,1)

    return(my_roc_plot)

}

# Construct h2o RF model
construct_rf <- function(
        train, # data to use for RF training
        cluster_label, # label of the cluster for which to develop the model
        response_column = "cluster_name", # column name containing labels
        nfolds = 5, # number of folds for K-fold cross-validation
        seed = 1 # seed to ensure reproducibility
    ){

    # Import a sample binary outcome train/test set into H2O
    new_labels <- ifelse(
        train[,response_column]==cluster_label,
        cluster_label,
        paste0("non",cluster_label)
    )

    train[,response_column] <- new_labels

    train <- as.h2o(train)

    # Identify predictors and response
    y <- response_column
    x <- setdiff(
        names(train),
        y
    )

    # For binary classification, response should be a factor
    train[, y] <- as.factor(train[, y])

    # Train & Cross-validate a RF
    my_rf <- h2o.randomForest(
        x = x,
        y = y,
        training_frame = train,
        ntrees = 50,
        nfolds = nfolds,
        keep_cross_validation_predictions = TRUE,
        seed = seed
    )

    return(my_rf)

}

# Plot votes confidence for a particular group
plot_votes <- function(
        class_to_plot # Which class(es) should be plotted
    ){

      votes_plot_data <- predictions %>%
          filter(
            predict %in% c(class_to_plot)
          )

      votes_plot <- votes_plot_data %>%
          ggplot(
            aes(
              x = as.character(pid),
              y = Confidence,
              color = Class,
              shape = Confidence <= 0.8
            )
          ) +
          ggbeeswarm::geom_quasirandom(size = 2) +
          GAMBLR::theme_Morons(
            my_legend_position = "right",
            my_legend_direction = "vertical",
            base_size = 22
          ) +
          geom_hline(
            yintercept = 0.8,
            linetype = "dashed",
            color = get_gambl_colours("lymphgenerator")[class_to_plot][1]
          ) +
          scale_color_manual(
            values = c(
              GAMBLR::get_gambl_colours("lymphgenerator"),
              GAMBLR::get_gambl_colours("lymphgen")
            )
          ) +
          ylab("Class assignment confidence") +
          ylim(c(0, 1)) +
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            )
          )

      return(votes_plot)
}


# initialize H2O
localH2O <- h2o.init(
    ip="localhost",
    port = 54321,
    startH2O = TRUE,
    nthreads=-1
)

# Prepare training (90%) and testing (10%) datasets
dataregression <- data_dlbcl %>%
    as.data.frame %>%
    t %>%
    as.data.frame %>%
    rownames_to_column("Tumor_Sample_Barcode") %>%
    left_join(
        .,
        res_meta %>%
          select(Tumor_Sample_Barcode, cluster_name),
        by = "Tumor_Sample_Barcode"
    ) %>%
    column_to_rownames(., "Tumor_Sample_Barcode")

set.seed(123)
training.samples <- dataregression$cluster_name %>%
    createDataPartition(
        p = 0.9,
        list = FALSE
    )
train.data <- dataregression[training.samples, ]
test.data <- dataregression[-training.samples, ]

# Number of cross-validation folds (to generate level-one data for stacking)
nfolds <- 5

# build RF models
rf_list <- lapply(
    unique(test.data$cluster_name),
    construct_rf,
    train = train.data
)

names(rf_list) <- unique(test.data$cluster_name)

# Now build the Ensembl model
# Import a sample binary outcome train/test set into H2O
train <- as.h2o(train.data)
test <- as.h2o(test.data)

# Identify predictors and response
y <- "cluster_name"
x <- setdiff(
    names(train),
    y
)

# For binary classification, response should be a factor
train[, y] <- as.factor(train[, y])
test[, y] <- as.factor(test[, y])

# Train a stacked ensemble using the RFs generated above
ensemble <- h2o.stackedEnsemble(
    x = x,
    y = y,
    training_frame = train,
    base_models = rf_list,
    seed = 1
)

h2o.performance(ensemble)

predictions <- h2o.predict(
        ensemble,
        newdata = test
    ) %>%
    as.data.frame %>%
    `rownames<-`(rownames(test.data)) %>%
    base::merge(
        .,
        test.data %>% select(cluster_name),
        by="row.names"
    ) %>%
    mutate(data_type="test") %>%
    rbind(
        .,
        (h2o.predict(ensemble, newdata = train) %>%
          as.data.frame %>%
          `rownames<-`(rownames(train.data)) %>%
          base::merge(
            .,
            train.data %>% select(cluster_name),
            by="row.names") %>%
            mutate(data_type="train")
        )
    ) %>%
    tidyr::pivot_longer(
        !c(Row.names, predict, cluster_name, data_type),
        names_to = "Class",
        values_to = "Confidence") %>%
    rename("pid" = "Row.names")

# Explore votes
plot_votes(c("aSEL", "aSCI"))
plot_votes(c("EGB", "ETB"))
plot_votes("MCaP")
plot_votes("MP3")

# which samples will have labels set to higher level due to lower confidence
# of their subgroup classification?
corrected_labels <- list()

corrected_labels$ST2 <- predictions %>%
    filter(
        predict %in% c("aSEL", "aSCI") |
        cluster_name %in% c("aSEL", "aSCI")
    ) %>%
    as.data.frame %>%
    group_by(pid) %>%
    mutate(
        first_vote = Confidence[which.max(Confidence)],
        second_vote = Confidence[which(rank(Confidence)==6)],
        second_class = Class[which(Confidence==second_vote)]
    ) %>%
    filter(first_vote <= 0.8) %>%
    distinct(pid) %>%
    pull(pid)

corrected_labels$EZB <- predictions %>%
  filter(
      predict %in% c("EGB", "ETB") |
      cluster_name %in% c("EGB", "ETB")
  ) %>%
  as.data.frame %>%
  group_by(pid) %>%
  mutate(
      first_vote = Confidence[which.max(Confidence)],
      second_vote = Confidence[which(rank(Confidence)==6)],
      second_class = Class[which(Confidence==second_vote)]
  ) %>%
  filter(first_vote <= 0.8) %>%
  distinct(pid) %>%
  pull(pid)

# Where would the corrected samples be on the phylo tree
# relative to the rest

d_ezb <- my_nj$tip.label %>%
  as.data.frame() %>%
  `colnames<-`("Tumor_Sample_Barcode") %>%
  left_join(
      .,
      umap_df %>%
        select(
          Tumor_Sample_Barcode,
          cluster_name
        )
      ) %>%
  mutate(
      label2 = ifelse(
        Tumor_Sample_Barcode %in% corrected_labels$EZB,
        Tumor_Sample_Barcode,
        ""
      )
  )

tree_ezb <- left_join(
        as.tibble(my_nj),
        d_ezb,
        by = c("label" = "Tumor_Sample_Barcode")
    ) %>%
    treeio::as.treedata()

phylo_plot_ezb <- ggtree::ggtree(
        tree_ezb,
        branch.length = "none",
        aes(
          color = cluster_name
        )
    ) +
    scale_color_manual(
        values = GAMBLR::get_gambl_colours("lymphgenerator"),
        labels = paste(
          "<span style='color:",
          GAMBLR::get_gambl_colours("lymphgenerator"),
          "'>",
          names(GAMBLR::get_gambl_colours("lymphgenerator")),
          "</span>",
          sep = ""
        )
    ) +
    labs(color = "Genetic subgroup") +
    ggtree::theme_tree2() +
    ggrepel::geom_label_repel(
        aes(label = label2),
        max.overlaps = 20,
        hjust = 0,
        size = 5,
        segment.size = 1,
        nudge_x = 5,
        direction = "y"
    ) +
    guides(
        color = guide_legend(override.aes = aes(label = ""))
    ) +
    theme(legend.text = ggtext::element_markdown())

phylo_plot_ezb

d_st2 <- my_nj$tip.label %>%
    as.data.frame() %>%
    `colnames<-`("Tumor_Sample_Barcode") %>%
    left_join(
        .,
        umap_df %>%
          select(
            Tumor_Sample_Barcode,
            cluster_name
          )
    ) %>%
    mutate(
        label2 = ifelse(
          Tumor_Sample_Barcode %in% corrected_labels$ST2,
          Tumor_Sample_Barcode,
          ""
        )
    )

tree_st2 <- left_join(
        as.tibble(my_nj),
        d_st2,
        by = c("label" = "Tumor_Sample_Barcode")
    )%>%
        treeio::as.treedata()

phylo_plot_st2 <- ggtree::ggtree(
        tree_st2,
        branch.length = "none",
        aes(
          color = cluster_name
        )
    ) +
    scale_color_manual(
        values = GAMBLR::get_gambl_colours("lymphgenerator"),
        labels = paste(
          "<span style='color:",
          GAMBLR::get_gambl_colours("lymphgenerator"),
          "'>",
          names(GAMBLR::get_gambl_colours("lymphgenerator")),
          "</span>",
          sep = ""
        )
    ) +
    labs(color = "Genetic subgroup") +
    ggtree::theme_tree2() +
    ggrepel::geom_label_repel(
        aes(label = label2),
        max.overlaps = 20,
        hjust = 0,
        size = 5,
        segment.size = 1,
        nudge_x = 12,
        direction = "y"
    ) +
    guides(
        color = guide_legend(override.aes = aes(label = ""))
    ) +
    theme(legend.text = ggtext::element_markdown())

phylo_plot_st2

###########################################################
####### Define the classification functionality ###########
###########################################################

# define final prediction function with high-level classes
predict_dlbcl <- function(
        h2o_model,
        list_of_base_models,
        new_data,
        high_level_group_cutoff = 0.8,
        unclassifiable_cutoff = 0.5,
        return_all_votes = TRUE
    ){

    # make a prediction
    prediction <- h2o.predict(
            h2o_model,
            newdata = as.h2o(new_data)
        ) %>%
        as.data.frame %>%
        `rownames<-`(
            rownames(new_data)
        )

    # tidy prediction
    prediction <- prediction %>%
        rownames_to_column("pid") %>%
        pivot_longer(
            !c(pid, predict),
            names_to = "Class",
            values_to = "Confidence"
        ) %>%
        group_by(pid) %>%
        arrange(
            desc(Confidence)
        ) %>%
        ungroup

    prediction <- prediction %>%
        group_by(pid) %>%
        mutate(
            predict = case_when(
                ((predict %in% c(
                    "aSEL",
                    "aSCI")
                  ) & Confidence <= high_level_group_cutoff) ~ "ST2",
                ((predict %in% c(
                    "EGB",
                    "ETB")
                  ) & Confidence <= high_level_group_cutoff) ~ "EZB",
                TRUE ~ Class[which.max(Confidence)])
        )

      # Review predictions and assign Other label when sample is ambiguous
      # First get samples in question (Confidence of vote < min cutoff)
      samples_in_question <- prediction %>%
          group_by(pid) %>%
          filter(
              Confidence == max(Confidence)
          ) %>%
          ungroup %>%
          filter(
              Confidence < unclassifiable_cutoff
          )

      # Now count how many featurees of the lead class they have
      # if there is only 1/top-10 for that class, sample is unclassifiable

      for(i in (1:nrow(samples_in_question))){

          # the sample id
          this_sample <- (samples_in_question[i, ] %>% pull(pid))

          # what sort of features this sample has?
          this_sample_features <- new_data[this_sample, ] %>%
              t %>%
              as.data.frame %>%
              filter(
                  !!as.name(this_sample) == 1
              )

          # which class this sample is predicted to be?
          this_prediction <- (samples_in_question[i, ] %>% pull(Class))

          # what are the features important for that class?
          this_prediction_important_feat <-
            pairwise_comparisons_all[[this_prediction]][["fisher"]] %>%
            filter(estimate > 1, q.value<0.1) %>%
            filter(gene %in%
                  (as.data.frame(
                    h2o.varimp(
                      list_of_base_models[[this_prediction]]
                    )
                  ) %>%
                  arrange(
                     desc(relative_importance)
                  ) %>%
                  pull(variable)
                  )
            ) %>%
            arrange(desc(estimate)) %>%
            head(10) %>%
            pull(gene)

          # what are the features from the class top-10 are in this sample?
          these_features <- rownames(
                  this_sample_features
              )[which(
                  rownames(
                    this_sample_features
                  ) %in% this_prediction_important_feat
                )]

          # sample is unclass if there is one or less features from that class
          if(length(these_features) < 2){
            prediction[prediction$pid == this_sample, "predict"] <- "Other"
          }
      }

      # convert predictions to factor so the related classes stay together
      order <- c(
          "MP3",
          "EZB", "EGB", "ETB",
          "ST2", "aSCI", "aSEL",
          "BNZ", "MCaP",
          "Other"
      )
      prediction$predict <- factor(
          prediction$predict,
          levels = order
      )

      # return only predicted label if desired
      if(!return_all_votes){
        prediction <- prediction %>%
          group_by(pid) %>%
          filter(
              Confidence==max(Confidence)
          ) %>%
          ungroup

        prediction <- prediction %>%
        arrange(match(
          pid,
          rownames(new_data)
        ))

        return(prediction)
      }

      prediction <- prediction %>%
        arrange(pid)

      return(prediction)

}

# Now make the predictions based on the model
predictions <- predict_dlbcl(
        ensemble,
        dataregression,
        list_of_base_models = rf_list,
        return_all_votes = FALSE
    ) %>%
    rename("Tumor_Sample_Barcode" = "pid") %>%
    left_join(
        .,
        res_meta
    )

###########################################################
####### Optimize the ensembl model ########################
###########################################################

# If we take the accuracy of predictions on the full data set
# as the ground truth, which features will decrease the accuracy
# of the prediction? By iterating through all features and leaving
# one of them out, let's build the new model and compare it's
# performance to the ground truth.

# Storage containers
metrics <- data.frame(
        Feature = character(),
        MSE = numeric(),
        Logloss = numeric(),
        stringsAsFactors = FALSE
    )

predictions <- data.frame(
        pid = character(),
        predict = character(),
        Class = character(),
        Confidence = numeric(),
        Feature = character(),
        stringsAsFactors = FALSE
    )


# Analyze all features
for(this_feature in colnames(dataregression %>% select(-cluster_name))){

    message(
        paste0(
          "Processing:",
          this_feature
        )
    )

    this_train.data <- train.data %>%
        select(-{{ this_feature }})
    this_test.data <- test.data %>%
        select(-{{ this_feature }})

    # build RF models
    this_rf_list <- list()

    for(group in unique(this_train.data$cluster_name)){
        print(group)
        this_rf <- construct_rf(
            train = this_train.data,
            cluster_label = group
        )
        this_rf_list[group] <- this_rf
    }

    rm(group)
    rm(this_rf)

    # Import a sample binary outcome train/test set into H2O
    this_train <- as.h2o(this_train.data)
    this_test <- as.h2o(this_test.data)

    # Identify predictors and response
    y <- "cluster_name"
    x <- setdiff(
            names(this_train),
            y
        )

    # For binary classification, response should be a factor
    this_train[, y] <- as.factor(this_train[, y])
    this_test[, y] <- as.factor(this_test[, y])

    # Train a stacked ensemble using the RF list above
    this_ensemble <- h2o.stackedEnsemble(
            x = x,
            y = y,
            training_frame = this_train,
            base_models = this_rf_list
        )

    these_metrics <- data.frame(
            Feature = this_feature,
            MSE = h2o.mse(this_ensemble),
            Logloss = h2o.logloss(this_ensemble)
        )

    metrics <- bind_rows(metrics, these_metrics)

    these_predictions <- predict_dlbcl(
            this_ensemble,
            bind_rows(
              this_train.data,
              this_test.data
            ),
            list_of_base_models = this_rf_list,
            return_all_votes = FALSE
        ) %>%
      mutate(Feature = this_feature)

    predictions <- bind_rows(
            predictions,
            these_predictions
        )

}

# Save the results
write_tsv(
        metrics,
        "results/loo_metrics.tsv"
    )

write_tsv(
        predictions,
        "results/loo_predictions.tsv"
    )

# What are the features that decrease the performance of the model?
kick_out <- c(metrics %>%
    filter(
        Logloss < h2o.logloss(ensemble),
        MSE < h2o.mse(ensemble)
    ) %>%
    arrange(Logloss) %>%
    pull(Feature),
    "FCGR2B_SV"
)

# Drop these features from the data
optimized_train.data <- train.data %>%
    select(-{{ kick_out }})
optimized_test.data <- test.data %>%
    select(-{{ kick_out }})


# build RF models
optimized_rf_list <- list()

optimized_rf_list <- lapply(
    unique(optimized_train.data$cluster_name),
    construct_rf,
    train = optimized_train.data
)

names(optimized_rf_list) <- unique(optimized_train.data$cluster_name)


# Import a sample binary outcome train/test set into H2O
optimized_train <- as.h2o(optimized_train.data)
optimized_test <- as.h2o(optimized_test.data)

# Identify predictors and response
x <- setdiff(
        names(optimized_train),
        y
    )

# For binary classification, response should be a factor
optimized_train[, y] <- as.factor(optimized_train[, y])
optimized_test[, y] <- as.factor(optimized_test[, y])

# Train a stacked ensemble using the RF list above
optimized_ensemble <- h2o.stackedEnsemble(
        x = x,
        y = y,
        training_frame = optimized_train,
        base_models = optimized_rf_list,
        seed = 1
    )

# Let's compare the votes between the base and optimized models
compare_votes <- bind_rows(
    predict_dlbcl(
        ensemble,
        bind_rows(
          train.data,
          test.data
        ),
        list_of_base_models = rf_list,
        return_all_votes = FALSE
    ) %>%
    select(pid, Confidence, predict) %>%
    mutate(method = "Production"),
    predict_dlbcl(
        optimized_ensemble,
        bind_rows(
          optimized_train.data,
          optimized_test.data
        ),
        list_of_base_models = optimized_rf_list,
        return_all_votes = FALSE
    ) %>%
    select(pid, Confidence, predict) %>%
    mutate(method = "Optimized"),
)

compare_votes %>%
    ggplot(
      aes(
        x = method,
        y = Confidence,
        color = predict
      )
    ) +
    geom_boxplot() +
    ggbeeswarm::geom_quasirandom() +
    facet_grid(
      rows = vars(predict), scales = "free_x"
    ) +
    theme_Morons() +
    ylab(
      "Vote confidence"
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      strip.text.y = element_text(angle = 0)
    ) +
    scale_color_manual(
      values = get_gambl_colours()
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = c(0, 0.01),
      limits = c(0, 1)
    )

ggpubr::ggpaired(
      compare_votes,
      x = "method",
      y = "Confidence",
      color = "predict",
      line.color = "predict",
      line.size = 0.2
    ) +
    facet_grid(
      cols = vars(predict),
      scales = "free_x"
    ) +
    theme_Morons() +
    scale_y_continuous(
      expand = c(0, 0.01),
      limits = c(0, 1)
    ) +
    ylab(
      "Vote confidence"
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle=90)
    ) +
    scale_color_manual(
      values = get_gambl_colours()
    )
