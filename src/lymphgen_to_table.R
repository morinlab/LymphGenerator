library(tidyverse)
library(tabulizer)

# locate areas only
f <- locate_areas("lymphgen.pdf")
f

# Extracting data
table <- extract_tables(
		"lymphgen.pdf",
		output = "data.frame",
		pages = 1, 
		area = list(
			c(10, 28, 551, 180)
		),
		guess = FALSE
	) %>% as.data.frame

# Minor wrangling
table <- table %>% # few columns are returned together
	separate(
		X.2,
		into = c("BN2", "N1", "EZB", "ST2"),
		sep = " "
	) %>%
	mutate( # because there are spaces in the column 2, create separator
		log = gsub(
			" -",
			"sep-",
			log
		)
	) %>% # and then split this column into 2
	separate(
		log,
		into = c("feature", "log_pval"),
		sep = "sep"
	)  %>% # rename the rest columns
	rename(
		"gene" = "X",
		"MCD" = "X.1",
		"A53" = "X.3"
	)

table <- table[-c(1,2),]
table <- table %>%
	mutate_at(
		c("log_pval", "MCD", "BN2", "N1", "EZB", "ST2", "A53"),
		as.numeric
	)
table$enriched_in <- c(
	rep("MCD", 38),
	rep("BN2", 21),
	rep("N1", 8),
	rep("EZB", 26),
	rep("ST2", 28),
	rep("A53", 46)
)

write_tsv(
	table,
	"data/lymphgen_classes.tsv"
)
