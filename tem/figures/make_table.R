

df <- tibble(file=list.files("/scratch/rheinnec/viktoria_figure/input/", pattern="correctionblend", recursive=T, full.names=T)) %>%
  mutate(
    label=c("A", "B", "C", "D"),
    start_x=c(5100, 2500, 2000, 1120),
    start_y=c(5100, 3000, 6000, 2000),
    width=c(1000, 2000, 3000, 4000)
  )


write_tsv(df, file="/scratch/rheinnec/viktoria_figure/df_in.tsv")



rel_files <- list.files("/g/schwab/tem_screen/organelles/images used/", pattern="*", recursive=T, full.names=F)



df_final <- tibble::tribble(
  ~label, ~start_x, ~start_y, ~width, ~id, 
  "BB",  2500,  3500,  3500,   "1_A5_Cut1_c016_116114586_TAL_20to200_20230628_TARA_epo_06_P2_correctionblend.png",
  "D",	2500,	3000,	4000,	"1_C1_Cut1_c004_116114875_BIL_20to200_20231007_TARA_epo_04_P2_correctionblend.png",
  "A",	1,	1,	3000,	"1_C1_Cut1_c005_116114875_BIL_20to200_20231007_TARA_epo_04_P2_correctionblend.png",
  "P",	4000,	4000,	5000,	"1_D6_Cut1_c021_116114875_BIL_20to200_20231007_TARA_epo_04_P1_correctionblend.png",
  "Y",	13000,	18000,	8000,	"245756_A5_Cut1_c001_116114425_TAL_10to40_20230617_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_A5_Cut1_c013_116114425_TAL_10to40_20230617_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_A5_Cut1_c015_116114425_TAL_10to40_20230617_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_A5_Cut1_c019_116114425_TAL_10to40_20230617_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_A5_Cut1_c020_116114425_TAL_10to40_20230617_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_A5_Cut1_c021_116114425_TAL_10to40_20230617_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_E5_Cut1_c001_117659905_BAR_20to200_20240321_AM_01_epo_03_P2_correctionblend.png",
  "A",	1,	1,	3000,	"245756_E5_Cut1_c003_117659905_BAR_20to200_20240321_AM_01_epo_03_P2_correctionblend.png",
  "A",	1,	1,	3000,	"245756_G1_Cut1_c003_117659905_BAR_20to200_20240321_AM_01_epo_03_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_H2_Cut1_c001_116115270_BAR_10to40_20240316_PM_01_epo_P3_correctionblend.png",
  "A",	1,	1,	3000,	"245756_H2_Cut1_c008_116115270_BAR_10to40_20240316_PM_01_epo_P3_correctionblend.png",
  "A",	1,	1,	3000,	"245756_J2_Cut2_c018_116114744_KRI_10to40_20230805_AM_01_epo_01_correctionblend.png",
  "A",	1,	1,	3000,	"245756_J2_Cut2_c019_116114744_KRI_10to40_20230805_AM_01_epo_01_correctionblend.png",
  "A",	1,	1,	3000,	"245756_K6_Cut1_c004_112518928_ATH_20to200_20240708_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_L1_Cut1_c008_116115038_NAP_10to40_20240416_AM_01_epo_03_P2_correctionblend.png",
  "A",	1,	1,	3000,	"245756_O2_Cut3_c005_116114983_POR_20to200_20231023_AM_01_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_P1_Cut1_c003_116114983_POR_20to200_20231023_AM_01_epo_01_P2_correctionblend.png",
  "A",	1,	1,	3000,	"245756_P1_Cut1_c004_116114983_POR_20to200_20231023_AM_01_epo_01_P2_correctionblend.png",
  "A",	1,	1,	3000,	"245756_Q4_Cut1_c026_116114998_NAP_10to40_20240415_AM_02_epo_01_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_R1_Cut1_c002_116115135_NAP_20to200_20240501_AM_01_epo_06_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_R1_Cut1_c003_116115135_NAP_20to200_20240501_AM_01_epo_06_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_S2_Cut2_c003_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_S2_Cut2_c004_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_S2_Cut2_c005_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c005_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c020_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c031_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c032_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c048_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c072_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c076_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png",
  "A",	1,	1,	3000,	"245756_T1_Cut1_c081_116114649_KRI_10to40_20230802_PM_01_epo_02_P1_correctionblend.png"
  
  
) %>%
  mutate(
    file=
      file.path(
        "/g/schwab/tem_screen/processed",
        str_remove(id, "_correctionblend.png"),
        str_replace(id, ".png$", ".mrc")
      )
    )



write_tsv(df_final[c(1,2,4,5),], file="/scratch/rheinnec/viktoria_figure/df_in.tsv")


write_tsv(df_final, file="/scratch/rheinnec/viktoria_figure/df_in_full.tsv")





library(tidyverse)
library(getopt)
library(googlesheets4)
library(googledrive)
library(cowplot)
#email = "marco.rheinnecker@embl.de"

json_key <- "/g/schwab/marco/repos/trec_biobank/tem/trec-tem-screen-e98a2e03f58b.json"
gs4_auth(path=json_key)
drive_auth(path = json_key)
trec_tem_googledoc <- "https://docs.google.com/spreadsheets/d/1VnX2JjlOJf7tkjw4DN6FIpiulz-qko1pQk6onMM5RwY/edit?gid=1366996799#gid=1366996799"



df <- read_sheet(trec_tem_googledoc, sheet="df_in_full", col_types="c") 


fig1 <- df  %>%
  filter(is.na(figure), !is.na(label))%>%
  mutate(label=factor(label, levels=c(LETTERS, paste0(LETTERS, LETTERS)))) %>%
  arrange(label)



#LETTERS


fig2 <- df %>% filter(figure==1) %>%
  mutate(label=factor(label, levels=c(LETTERS, paste0(LETTERS, LETTERS)))) %>%
  arrange(label)



write_tsv(fig1, file="/scratch/rheinnec/viktoria_figure/df_in_full_new_big.tsv")


write_tsv(fig2, file="/scratch/rheinnec/viktoria_figure/df_in_full_new.tsv")
















