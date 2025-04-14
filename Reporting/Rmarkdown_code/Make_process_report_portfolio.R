library(stringi)
library(lubridate) 
library(knitr) # to make html tables
library(kableExtra) # to make html tables
library(dplyr)
library(stringr)
library(rmarkdown)

print("loaded libraries")

#Document header for HTML report
document_number = "Not Cleared"
version_number = " - Draft Document"
doc_control <- paste("Document Control: ",document_number, version_number, sep="")
head_doc_control <- doc_control

#Get Date Info
time_now <- format(Sys.time(), "%H%M")
date_now <- Sys.Date()

## DATE INFORMATION
current_date <- paste("Date (YYYY-MM-DD): ",date_now,"
", sep="")
fullPath <- "/Users/davejacobson/Desktop/DataScience_Trainings/Portfolio/Reporting/Rmarkdown_code/"
epi_clus <- read.csv(paste(fullPath,"input_data/epi_clusters.csv", sep  =""), header = T)
lab_data <- read.csv(paste(fullPath,"input_data/lab_metadata.csv",sep  =""), header =T)
clus_data <- read.csv(paste(fullPath,"input_data/cluster_memberships.tsv",sep  =""), sep = "\t", header =T)
tgc_data <- read.csv(paste(fullPath,"input_data/tgc_mem.csv",sep  =""), sep = ",", header =T)

#remove columns that have empty header names
epi_clus <- epi_clus[,!grepl("X", colnames(epi_clus))]
lab_data <- lab_data[,!grepl("X", colnames(lab_data))]
clus_data <- clus_data[,!grepl("X", colnames(clus_data))]
tgc_data <- tgc_data[,!grepl("X", colnames(tgc_data))]

lab_plus_epi <- merge(lab_data, epi_clus, by.x = "Case_ID", by.y = "Case_ID", all.x = T)
withClus <- merge(lab_plus_epi, clus_data, by.x = "LSDB_Sequence_ID", by.y = "Seq_ID", all.x = T)
latestMaster <- merge(withClus, tgc_data, by.x = "LSDB_Sequence_ID", by.y = "LSDB_Sequence_ID", all.x = T)

#fixing some column names to match
current_clinical_cases <- latestMaster
names(current_clinical_cases)[names(current_clinical_cases) == "LSDB_Sequence_ID"] <- "LSDB_Seq_ID"
names(current_clinical_cases)[names(current_clinical_cases) == "LSDB_Specimen_ID"] <- "LSDB_Lab_ID"
names(current_clinical_cases)[names(current_clinical_cases) == "Date_Collected"] <- "Collection_date"
names(current_clinical_cases)[names(current_clinical_cases) == "Cluster"] <- "GeneticCluster"
current_clinical_cases$Date_changed_TGC[is.na(current_clinical_cases$Date_changed_TGC)] <- "4/7/2025"
# print(current_clinical_cases)
# quit()
#adding in residence column here for example purposes
current_clinical_cases$State_of_Residence <- current_clinical_cases$Submitter_State

#Sort by seq id so that tables are ordered by LSDB ID for each state
current_clinical_cases <- current_clinical_cases[order(current_clinical_cases$LSDB_Seq_ID),]

#Sequencing states won't have an LSDB Lab ID but they will have a Sequence ID. Populate Lab ID with Seq ID for these states
current_clinical_cases$LSDB_Lab_ID[is.na(current_clinical_cases$LSDB_Lab_ID)] <- current_clinical_cases$LSDB_Seq_ID
current_clinical_cases$LSDB_Lab_ID[current_clinical_cases$LSDB_Lab_ID == ""] <- current_clinical_cases$LSDB_Seq_ID
current_clinical_cases$LSDB_Lab_ID[current_clinical_cases$LSDB_Lab_ID == " "] <- current_clinical_cases$LSDB_Seq_ID

#Get a list of all the states that have submitted samples and have a resident with a sample. Combine those lists and start making reports
all_states_to_make_a_report_for1 <- unique(current_clinical_cases$Submitter_State)
all_states_to_make_a_report_for2 <- unique(current_clinical_cases$State_of_Residence)
all_states_to_make_a_report_for3 <- c(all_states_to_make_a_report_for1 ,all_states_to_make_a_report_for2)
all_states_to_make_a_report_for_final <- unique(all_states_to_make_a_report_for3)
#If there are reports to make for any state
if(length(all_states_to_make_a_report_for_final) > 0){ 

	#loop through each state
	for( k in all_states_to_make_a_report_for_final){
		if(k != ""){
			print("start state")
			print(k)
		

	current_state_code <- k

	## STATE INFORMATION for HTML report
	current_state_running <- paste0("State Partner: ",current_state_code)

	#Add underscore before State Abbreviation to make it easier to pull out the state 
	state_string <- paste0("_",k)

	CURRENT_clinical_cases_state_dt <- current_clinical_cases %>% filter(str_detect(Submitter_State, paste(k,"$", sep ="")))

	#Check status of positive and negative sequencing controls. Make a table for sequencing states, non-sequencing states table is left blank
	PosDF <- CURRENT_clinical_cases_state_dt[grepl("PS", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]
	NegDF <- CURRENT_clinical_cases_state_dt[grepl("NG", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]

	#Grab the row number if a positive control has failed
	pos_fail_vec <- c()
	if(length(rownames(PosDF))> 0){
		state_df_pos_final <- PosDF[,c(1,13,14)]
		for(a in 1:length(rownames(PosDF))){
			if(state_df_pos_final[a,2] == "FAIL"){
				pos_fail_vec <- c(pos_fail_vec, a)
			}
		}
	} else{
		state_df_pos_final <- data.frame(matrix(ncol = 3, nrow = 1))
		state_df_pos_final[1,] <- "N/A"
	}

	#For negatives, we have to switch the FAIL/PASS results. Since a "PASS" means it has a succesfully cluster genotype, we want this to be a 'FAIL' for negative controls
	#Also grab the row number of the negative controls that failed
	neg_fail_vec <- c()
	if(length(rownames(NegDF))> 0){
		state_df_neg_final <- NegDF[,c(1,13,14)]
		for(a in 1:length(rownames(NegDF))){
			if(state_df_neg_final[a,2] == "FAIL"){
				print("change to pass")
				state_df_neg_final[a,2] <- "PASS"
		} else if(state_df_neg_final[a,2] == "PASS"){
			print("change to fail")
			state_df_neg_final[a,2] <- "FAIL"
			neg_fail_vec <- c(neg_fail_vec, a)
		}
		}
	} else{
		state_df_neg_final <- data.frame(matrix(ncol = 3, nrow = 1))
		state_df_neg_final[1,] <- "N/A"
	}


	rownames(state_df_neg_final) <- NULL
	rownames(state_df_pos_final) <- NULL
	colnames(state_df_pos_final) <- c("Seq_ID", "Status of positive control", "Reason for specimen failure (if applicable)")
	colnames(state_df_neg_final) <- c("Seq_ID", "Status of negative control", "Haplotypes not detected in negative control (if applicable)")

	#Make the pos/neg tables for reports
	pos_controls <- knitr::kable(state_df_pos_final, "html", padding = 50, line_sep = 2, table.attr = "style='width:50%;'", align = rep("c", 3)) %>%  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "float_left", font_size = 12)
	pos_controls <- add_header_above(pos_controls, c("Results obtained for Positive controls (PS)" = 3), italic = FALSE, font_size = 14, align = "center", line = FALSE, background = "#9A9A9A", color = "#F5F5F5", extra_css = "line-height: 20pt")
	pos_controls <- footnote(pos_controls, general_title= "Note", title_format = "bold", general = "Positive controls are considered PASS if they pass the Cyclospora genotyping inclusion criteria. Positive controls are considered FAIL if they do not pass the Cyclospora genotyping inclusion criteria.")

	neg_controls <- knitr::kable(state_df_neg_final, "html", padding = 50, line_sep = 2, table.attr = "style='width:50%;'", align = rep("c", 3)) %>%  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F, position = "left", font_size = 12)
	neg_controls <- add_header_above(neg_controls, c("Results obtained for Negative controls (NG)" = 3), italic = FALSE, font_size = 14, align = "center", line = FALSE, background = "#9A9A9A", color = "#F5F5F5", extra_css = "line-height: 20pt")
	neg_controls <- footnote(neg_controls, general_title= "Note", title_format = "bold", general = "Negative controls are considered PASS if the do not pass the Cyclospora genotyping inclusion criteria. Negative controls are considered FAIL if they pass the Cyclospora genotyping inclusion criteria.")

	#Add color for pos/neg rows that are fails
	if(length(neg_fail_vec) > 0){
		neg_controls <- row_spec(neg_controls,neg_fail_vec, background = "#ff8080", align = "c")
	}

	if(length(pos_fail_vec) > 0){
		pos_controls <- row_spec(pos_controls,pos_fail_vec, background = "#FFEC8B", align = "c")
	}


	#Remove controls from the dataset for making tables for clinical samples
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[!grepl("NG", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[!grepl("PS", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[!grepl("VD", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]

	#Make sure the Lab IDs are filled in with Seq IDs (should be done before the loop but this is a double check)
	CURRENT_clinical_cases_state_dt$LSDB_Lab_ID[is.na(CURRENT_clinical_cases_state_dt$LSDB_Lab_ID)] <- CURRENT_clinical_cases_state_dt$LSDB_Seq_ID
	CURRENT_clinical_cases_state_dt$LSDB_Lab_ID[CURRENT_clinical_cases_state_dt$LSDB_Lab_ID == ""] <- CURRENT_clinical_cases_state_dt$LSDB_Seq_ID
	CURRENT_clinical_cases_state_dt$LSDB_Lab_ID[CURRENT_clinical_cases_state_dt$LSDB_Lab_ID == " "] <- CURRENT_clinical_cases_state_dt$LSDB_Seq_ID


	#extract columns that will be used for different purposes
	#This small table is used for the TGC updates, so don't need everything
	print(CURRENT_clinical_cases_state_dt)
	# quit()
	CURRENT_clinical_cases_state_dt_smallTable <- CURRENT_clinical_cases_state_dt[,c(1,2,4,5,6,9,10,13,14,15,16)]
	# CURRENT_clinical_cases_state_dt_smallTable <- CURRENT_clinical_cases_state_dt[,c(1:16)]

	#this next table doesn't include the previous TGCs and it is used for making the table with all specimens
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[,c(1:16)]

	#this subset table is used for making sure we get all specimens pertaining to the state (regardless of who submitted them)
	current_clinical_cases_subset <- current_clinical_cases[,c(1:16)]
	CURRENT_clinical_cases_submitted_by_otherState <- current_clinical_cases_subset %>% filter(str_detect(State_of_Residence,paste(k, "$", sep ="")))
	CURRENT_clinical_cases_state_dt <- rbind(CURRENT_clinical_cases_state_dt,CURRENT_clinical_cases_submitted_by_otherState)

	#Remove any duplicated rows
	#CURRENT_clinical_cases_state_dt <- distinct(CURRENT_clinical_cases_state_dt, LSDB_Lab_ID, .keep_all =T)
	CURRENT_clinical_cases_state_dt  <- CURRENT_clinical_cases_state_dt[!duplicated(CURRENT_clinical_cases_state_dt),]

	#Drop controls from combined table of samples submitted from the state and samples from residents of the states
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[!grepl("NG", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[!grepl("PS", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[!grepl("VD", CURRENT_clinical_cases_state_dt$LSDB_Seq_ID),]

	#fill in 'not provided' or 'N/A' for missing data, where applicable
	CURRENT_clinical_cases_state_dt$Case_ID[CURRENT_clinical_cases_state_dt$Case_ID=="None"] <- "not provided"
	CURRENT_clinical_cases_state_dt$Case_ID[CURRENT_clinical_cases_state_dt$Case_ID==""] <- "not provided"
	CURRENT_clinical_cases_state_dt$Case_ID[is.na(CURRENT_clinical_cases_state_dt$Case_ID)] <- "not provided"
	CURRENT_clinical_cases_state_dt$Collection_date[CURRENT_clinical_cases_state_dt$Collection_date=="None"] <- "not provided"
	CURRENT_clinical_cases_state_dt$Collection_date[CURRENT_clinical_cases_state_dt$Collection_date==""] <- "not provided"
	CURRENT_clinical_cases_state_dt$Collection_date[is.na(CURRENT_clinical_cases_state_dt$Collection_date)] <- "not provided"

	CURRENT_clinical_cases_state_dt$Reason_for_fail[CURRENT_clinical_cases_state_dt$Reason_for_fail=="None"] <- "N/A"
	CURRENT_clinical_cases_state_dt$Reason_for_fail[CURRENT_clinical_cases_state_dt$Reason_for_fail==""] <- "N/A"
	CURRENT_clinical_cases_state_dt$Reason_for_fail[is.na(CURRENT_clinical_cases_state_dt$Reason_for_fail)] <- "N/A"


	#Calculate how many specimens were submitted by state but resident of other state, etc
	submitted_total <- length(CURRENT_clinical_cases_state_dt$LSDB_Lab_ID)
	submitted_and_resident <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(State_of_Residence, paste(k,"$", sep =""))) %>% filter(str_detect(Submitter_State, paste(k,"$", sep =""))) %>% count(Submitter_State)
	submitted_not_resident <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Submitter_State, paste(k,"$", sep =""))) %>% filter(!State_of_Residence %in% paste(k,"$", sep ="")) %>% count(Submitter_State)
	submitted_byOtherState_but_resident <- CURRENT_clinical_cases_state_dt %>% filter(State_of_Residence %in% paste(k,"$", sep ="")) %>% filter(!Submitter_State %in% paste(k,"$", sep ="")) %>% count(State_of_Residence)

	#I think we only need the assignment or the if/else statements, not both. Can verify later
	submitted_and_resident_val <- 0
	submitted_not_resident_val <- 0

	submitted_and_resident_val <- submitted_and_resident$n
	submitted_not_resident_val <- submitted_not_resident$n

	if(length(rownames(submitted_and_resident)) >0 ){
		submitted_and_resident_val <- submitted_and_resident$n
	}else{
		submitted_and_resident_val <- 0
	}

	if(length(rownames(submitted_not_resident)) >0 ){
		submitted_not_resident_val <- submitted_not_resident$n
	}else{
		submitted_not_resident_val <- 0
	}

	if(length(rownames(submitted_byOtherState_but_resident)) >0 ){
		submitted_byOtherState_but_resident_val <- submitted_byOtherState_but_resident$n
	}else{
		submitted_byOtherState_but_resident_val <- 0
	}

	#Pass / Fail / Pending per state
	CURRENT_made_it_to_tree <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Status_of_specimen, "PASS"))
	CURRENT_NUMBER_success <- length(CURRENT_made_it_to_tree$LSDB_Seq_ID)

	CURRENT_failed_to_make_tree <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Status_of_specimen, "FAIL"))
	CURRENT_NUMBER_failed <- length(CURRENT_failed_to_make_tree $LSDB_Lab_ID)

	CURRENT_pending_status <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Status_of_specimen, "PENDING"))
	CURRENT_NUMBER_pending <- length(CURRENT_pending_status $LSDB_Lab_ID)

	submitted_for_this_state <- (CURRENT_NUMBER_success + CURRENT_NUMBER_failed + CURRENT_NUMBER_pending)

	#Metadata missing for any cases - use this info for bringing missing data to attention of states
	CURRENT_missing_Case_ID <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Case_ID, "not provided"))
	CURRENT_count_missing_Case_ID <- length(CURRENT_missing_Case_ID$CDC_Lab_ID)

	CURRENT_missing_Collection_date <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Collection_date, "not provided"))
	CURRENT_count_missing_Collection_date <- length(CURRENT_missing_Collection_date$CDC_Lab_ID)


	#Format so that all Not Applicable all look the same
	CURRENT_clinical_cases_state_dt$GeneticCluster <- gsub("NA", "N/A", CURRENT_clinical_cases_state_dt$GeneticCluster)

	#Check if most recent status change is less than 7 days, this means the specimen should be highlighted for a reason
	#Can either because it is a new specimen (added to report = status change), or TGC Change (added to report != status change)
	#This code can be buggy if the date is formatted incorrectly. May need to check this loop, print results, and fix samples in the working master
	# Check if status change equals date added to report - that means it is a new specimen
	#I put these into vectors that can later be used to specify which rows should be highlighted
	newSpecimenVec <- c()
	statusChangeVec <- c()
	noRemarkVec <- c()
# print(CURRENT_clinical_cases_state_dt)
	if(length(rownames(CURRENT_clinical_cases_state_dt))==0){
		for(i in 16){
			CURRENT_clinical_cases_state_dt[1,i] <- "N/A"
		}
	}else{
	for(a in 1:length(rownames(CURRENT_clinical_cases_state_dt))){
		addedToReport <- as.Date(CURRENT_clinical_cases_state_dt[a,4], "%m/%d/%Y")
		statusChange <- as.Date(CURRENT_clinical_cases_state_dt[a,15], "%m/%d/%Y")
		# print(addedToReport)
		# print(statusChange)
		timeDiff <- date_now - statusChange
		# print(CURRENT_clinical_cases_state_dt[a,])
		# print(statusChange)
		# print(timeDiff)
		if(timeDiff < 7){
			if(addedToReport != statusChange){
				#This will be only for samples that change TGCs (will have the TGC code output date an update was made)
				#This should status changes that are not due to being a new specimen
				#16 should be a new column where we can track what we want to highlight
				CURRENT_clinical_cases_state_dt[a,16] <- "StatusChange"
				statusChangeVec <- c(statusChangeVec,a)
			}
			else{ 
				CURRENT_clinical_cases_state_dt[a,16] <- "NewSpecimen"
				newSpecimenVec <- c(newSpecimenVec, a)
			}
		}
		else{
			#If status change is more than 7 days ago, we won't highlight it becuaes it is not noteworthy
			CURRENT_clinical_cases_state_dt[a,16] <- "NoRemark"
			noRemarkVec <- c(noRemarkVec,a)
		}

	}
	}
	
	colnames(CURRENT_clinical_cases_state_dt)[16] <- "reportHighligting"
	# print(CURRENT_clinical_cases_state_dt)

	#we'll come back to working with the status change specimens later on

	#Reorder columns so they fit a more logical flow of information
	
	# CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[,c(7, 8, 9, 2, 3,4,5,6, 13, 14, 10, 11, 12, 15, 1, 16)]
	CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[,c(4,4,15,2,6,5,7,7,9,8,13,10,1)]

	#Drop the column that has the statuschange / new specimen information for specimens, because we are using row indices to color those rows
	# CURRENT_clinical_cases_state_dt <- CURRENT_clinical_cases_state_dt[,-16]
	CURRENT_clinical_cases_state_dt_forKable <- CURRENT_clinical_cases_state_dt

	#Drop genetic cluster and lsdb lab ID, this can change based on what we decide we want to include in reports
	# CURRENT_clinical_cases_state_dt_forKable <- CURRENT_clinical_cases_state_dt_forKable[,-c(11,14)]


	#Rownames should be null for making kable table
	rownames(CURRENT_clinical_cases_state_dt_forKable) <- NULL
	#Rename columns to be a little bit shorter in an effort to keep the table as condensed as possible
	colnames(CURRENT_clinical_cases_state_dt_forKable) <- c("Received", "First_Reported", "Updated","Case_ID", "State_Lab_ID", "Collected", "Submitter", "Residence", "Status","Reason_for_fail^1^", "TGC_Code^2^", "Cluster","LSDB_Seq_ID")
	# print(CURRENT_clinical_cases_state_dt_forKable)
	#Make the kable for the HTML report
	kable_out <- knitr::kable(CURRENT_clinical_cases_state_dt_forKable, "html", padding = 10, line_sep = 2, align = rep("c", 13)) %>% 
		kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
			column_spec(c(3,8,11,12), border_right = "2px solid gray")

	#Highlight rows based on if specimens have changed TGCs or are new
	if(length(newSpecimenVec > 0)){
		kable_out <- row_spec(kable_out,newSpecimenVec, background = "#FFEC8B", align = "c")
		newSpecimenTable <- CURRENT_clinical_cases_state_dt_forKable[newSpecimenVec,] 
	}else{
		newSpecimenTable <- data.frame(matrix(ncol = 13, nrow = 1))
		colnames(newSpecimenTable) <- colnames(CURRENT_clinical_cases_state_dt_forKable)
		#This means that there aren't specimens that changed tgcs. Still make the table, but fill with N/As
		newSpecimenTable[1,] <- "N/A"
	}
	if(length(statusChangeVec > 0)){
		kable_out <- row_spec(kable_out,statusChangeVec, background = "#F88D76", align = "c")
		#If there are specimens that change TGCs, then make a table for them
		statusChangeTable <- CURRENT_clinical_cases_state_dt_smallTable[statusChangeVec,]
	} else{
		statusChangeTable <- data.frame(matrix(ncol = 11, nrow = 1))
		colnames(statusChangeTable) <- colnames(CURRENT_clinical_cases_state_dt_smallTable)
		#This means that there aren't specimens that changed tgcs. Still make the table, but fill with N/As
		statusChangeTable[1,] <- "N/A"
	}

	#Start adding headers to the main table
	kable_out <- add_header_above(kable_out, c("CDC Processing Timeline"= 3, "Specimen Metadata" = 5, "Genotyping Results" = 3, "Epi" = 1, "CDC ID" = 1), italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", background = "#FAEBD7", color = "black", extra_css = "line-height: 20pt;")
	kable_out <- add_header_above(kable_out, c("Specimens Genotyped"= 13), italic = FALSE, font_size = 16, align = "center", line = FALSE, background = "#9A9A9A", color = "#F5F5F5", extra_css = "line-height: 20pt;")

	#Collect specimens with missing metadata.
	those_with_missing_dates <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Collection_date, "not provided"))
	those_with_missing_case_ids <- CURRENT_clinical_cases_state_dt %>% filter(str_detect(Case_ID, "not provided"))
			
	list_with_missing_dates <- those_with_missing_dates$State_Lab_ID 
	list_with_missing_case_ids <- those_with_missing_case_ids$State_Lab_ID    	

	#They need to be put into a string to add the top of the table
	specimens_needing_dates <- toString(list_with_missing_dates)
	specimens_needing_case <- toString(list_with_missing_case_ids)

	#These strings need to have a new line bewteen entries to look nice for the reports
	specimens_needing_dates <- str_replace_all(specimens_needing_dates, ",", "
	")
	specimens_needing_case <- str_replace_all(specimens_needing_case, ",", "
	")

	missing_dates_string <- paste("Please provide CDC with collection dates for these specimens: 

	", specimens_needing_dates, sep="")
	missing_case_ids_string <- paste("Please provide CDC with Case IDs for these specimens: 

	", specimens_needing_case, sep="")

	#Make the missing data strings formatted to fit the table
	missing_data <- c("","","",missing_dates_string,"","","","", missing_case_ids_string,"","","","")

	#Add the missing metadata info to the top of the main table
	kable_out <- add_header_above(kable_out, c(" "= 13), italic = FALSE, font_size = 8, align = "justify", line = FALSE, color = "white", background = "#F5F5F5") %>%
	add_header_above(missing_data, italic = FALSE, font_size = 12, align = "c", line = FALSE, color = "red", extra_css = "line-height: 16pt; border-bottom: 1px dotted black; vertical-align:top;") #can also align to 'top' and 'bottom'

	kable_out <- add_header_above(kable_out, c(" "= 13), italic = FALSE, font_size = 8, align = "justify", line = FALSE, color = "white", background = "#F5F5F5")

	kable_out <- add_header_above(kable_out, c("IMPORTANT: If there are specimens for which CDC requires  additional information, IDs for these  specimens will be listed in red text below. Please provide CDC with this information or these specimens cannot be processed any further. Forward the relevant information as a line list in excel or csv format, including the Case ID and/or Collection date for these specimens to the following email address: CyclosporaAMD@cdc.gov"= 13), italic = FALSE, font_size = 14, align = "justify", line = TRUE, color = "red", background = "white") 
	kable_out <- add_header_above(kable_out, c(" "= 13), italic = FALSE, font_size = 12, align = "justify", line = TRUE, color = "white", background = "#F5F5F5")

	kable_out <- add_header_above(kable_out, c("Additional Information Required"= 13), italic = FALSE, font_size = 16, align = "center", line = FALSE, color = "#F5F5F5", background = "#9A9A9A", extra_css = "line-height: 20pt")
	#make footnote
	kable_out <- footnote(kable_out, number=c("A specimen may be sequenced multiple times if the first sequencing run failed. If the second (or third) sequencing of a specimen is successful it will be listed as PASS but the reason for fail column will still print the failed markers from the previous sequencing run. This is simply for reference and should not be a cause for concern.", "This report contains limited information on the temporal genetic clusters detected using these genotyped specimens. For more information please consult your state epidemologists, Cyclospora@cdc.gov, and/or CyclosporaAMD@cdc.gov"))


	#Format the kable column widths
	if(length(CURRENT_clinical_cases_state_dt$LSDB_Lab_ID) > 0) {
		
	kable_out <- column_spec(kable_out, 1, width = "4cm")
	kable_out <- column_spec(kable_out, 2, width = "4cm")
	kable_out <- column_spec(kable_out, 3, width = "4cm")

	kable_out <- column_spec(kable_out, 4, width = "4cm")
	kable_out <- column_spec(kable_out, 5, width = "4cm")
	kable_out <- column_spec(kable_out, 6, width = "4cm")

	kable_out <- column_spec(kable_out, 7, width = "4cm")
	kable_out <- column_spec(kable_out, 8, width = "4cm")
	kable_out <- column_spec(kable_out, 9, width = "4cm")
	kable_out <- column_spec(kable_out, 10, width = "4cm")

	kable_out <- column_spec(kable_out, 11, width = "4cm")
	kable_out <- column_spec(kable_out, 12, width = "4cm")
	kable_out <- column_spec(kable_out, 13, width = "4cm")
	# kable_out <- column_spec(kable_out, 14, width = "4cm")
	# kable_out <- column_spec(kable_out, 15, width = "4cm")
	# kable_out <- column_spec(kable_out, 16, width = "4cm")

	}


	#new data frames for the results summary. Use number calculated earlier on in the script
	summary_data_frame = data.frame(matrix("", ncol = 2, nrow = 4)) 

	summary_data_frame[1,1] <- "Total number of specimens submitted for genotyping (this includes all specimens submitted by the state lab plus any specimens submitted for local case patients by out-of-state labs):  "
	summary_data_frame[1,2] <- submitted_for_this_state

	summary_data_frame[2,1] <- "Number of specimens submitted by the state lab for local case patients:  "
	summary_data_frame[2,2] <- submitted_and_resident_val

	summary_data_frame[3,1] <- "Number of specimens submitted by the state lab for out-of-state case patients:  "
	summary_data_frame[3,2] <- submitted_not_resident_val

	summary_data_frame[4,1] <- "Number of specimens submitted by out-of-state lab(s) for local case patients:  "
	summary_data_frame[4,2] <- submitted_byOtherState_but_resident_val


	#slightly different summary table
	firstTable_df <- data.frame(matrix("", ncol = 2, nrow = 4)) 

	firstTable_df[1,1] <- "Total number of specimens submitted for genotyping (this includes all specimens submitted by the state lab plus any specimens submitted for local case patients by out-of-state labs):  "
	firstTable_df[1,2] <- submitted_for_this_state

	firstTable_df[2,1] <- "Number of successfully genotyped specimens (i.e., specimens passing our inclusion criteria):  "
	firstTable_df[2,2] <- CURRENT_NUMBER_success

	firstTable_df[3,1] <- "Number of specimens that failed to genotype (see table below for reason each specimen failed genotyping):  "
	firstTable_df[3,2] <- CURRENT_NUMBER_failed

	firstTable_df[4,1] <- "Number of specimens with status pending (e.g., awaiting genotyping):  "
	firstTable_df[4,2] <- CURRENT_NUMBER_pending

	rownames(summary_data_frame) <- NULL
	colnames(summary_data_frame) <- NULL
	rownames(firstTable_df) <- NULL
	colnames(firstTable_df) <- NULL


	#Make the summary tables
	SUMMARY <- knitr::kable(summary_data_frame, "html", padding = 10, line_sep = 2, table.attr = "style='width:100%;'", align = c("l", "r")) %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12)
	SUMMARY <- add_header_above(SUMMARY, c("Specimen Processing Information (only for current cases - excludes controls and validation specimens)
	"= 2), italic = FALSE, font_size = 20, align = "center", line = FALSE, background = "#9A9A9A", color = "#F5F5F5", extra_css = "line-height: 20pt;")

	firstTable <- knitr::kable(firstTable_df, "html", padding = 10, line_sep = 2, table.attr = "style='width:100%;'", align = c("l", "r")) %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12)
	firstTable <- add_header_above(firstTable, c("Specimen Processing Status for All Submitted Specimens"= 2), italic = FALSE, font_size = 20, align = "center", line = FALSE, background = "#9A9A9A", color = "#F5F5F5", extra_css = "line-height: 20pt;")

	#Add a footnote
	firstTable <- footnote(firstTable, general_title = "Note", title_format = "bold", general = "If you feel that the contents of the above tables do not reflect the correct number of submitted/processed specimens please notify the CDC *Cyclospora* genotyping laboratory at CyclosporaAMD@cdc.gov")

	rownames(statusChangeTable) <- NULL
	#Reorganize the tgc changed table to follow similar flow of main table
	print(statusChangeTable)
	# statusChangeTable <- statusChangeTable[,c(7,2,3,4,5,6,10,8,11,9,1)]
	statusChangeTable <- statusChangeTable[,c(10,2,5,4,11,11,9,8,6,7,1)]
	colnames(statusChangeTable) <- c("Updated", "Case_ID", "State_Lab_ID", "Collection_Date", "Submitter", "Residence", "Prev. TGC", "Current TGC", "Status", "Cluster", "LSDB_Seq_ID")
	statusChangeTable_Out <- knitr::kable(statusChangeTable, "html", padding = 10, line_sep = 2, align = rep("c", 11)) %>% kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12)

	#Width doesn't seem to be changed with this code, so probably could be deleted. But do keep the column outlining
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, c(1,6,9,10), border_right = "2px solid gray")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 1, width = "4cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 2, width = "5cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 3, width = "5cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 4, width = "5cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 5, width = "4cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 6, width = "4cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 7, width = "5cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 8, width = "5cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 9, width = "8cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 10, width = "8cm")
	statusChangeTable_Out <- column_spec(statusChangeTable_Out, 11, width = "6cm")

	#Format the table
	statusChangeTable_Out <- add_header_above(statusChangeTable_Out, c("Timeline"= 1, "Specimen Metadata" = 5, "TGC Data" = 3, "Epi Data" = 1, "CDC ID" = 1), italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", background = "#FAEBD7", color = "black", extra_css = "line-height: 20pt;")
	statusChangeTable_Out <- add_header_above(statusChangeTable_Out, c("Specimens with a TGC Status Change over the Past Week"=11), italic = FALSE, font_size = 20, align = "center", line = FALSE, background = "#9A9A9A", color = "#F5F5F5", extra_css = "line-height: 20pt;")

	rownames(newSpecimenTable) <- NULL
	colnames(newSpecimenTable) <- c("Received", "First_Reported", "Updated","Case_ID", "State_Lab_ID", "Collected", "Submitter", "Residence", "Status","Reason_for_fail", "TGC_Code", "Cluster","LSDB_Seq_ID")
	newSpecimenTable <- newSpecimenTable[,-3]
	newSpecimen_kable_out <- knitr::kable(newSpecimenTable, "html", padding = 10, line_sep = 2, align = rep("c", 12)) %>% 
		kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = T, position = "center", font_size = 12) %>%
			column_spec(c(2,7,10,11), border_right = "2px solid gray")

	newSpecimen_kable_out <- add_header_above(newSpecimen_kable_out, c("CDC Processing Timeline"= 2, "Specimen Metadata" = 5, "Genotyping Results" = 3, "Epi" = 1, "CDC ID" = 1), italic = FALSE, font_size = 14, align = "center", line = T, border_left = T, border_right = "2px solid gray", background = "#FAEBD7", color = "black", extra_css = "line-height: 20pt;")
	newSpecimen_kable_out <- add_header_above(newSpecimen_kable_out, c("Specimens Genotyped"= 12), italic = FALSE, font_size = 16, align = "center", line = FALSE, background = "#9A9A9A", color = "#F5F5F5", extra_css = "line-height: 20pt;")

	#Make the table and name it / put it in the right spot
	rename_html <- paste0(fullPath,"output_reports/",date_now,"_",time_now,"_",current_state_code,"_cycloGenotyping_process_report.html")
	rmarkdown::render(paste(fullPath,"Make_process_report_portfolio.Rmd",sep =""), "html_document",  run_pandoc = TRUE)

	file.copy(paste(fullPath,"Make_process_report_portfolio.html", sep = ""), to= rename_html, overwrite = TRUE)
	file.remove(paste(fullPath,"Make_process_report_portfolio.html",sep = ""))

	}
	}

}








 

