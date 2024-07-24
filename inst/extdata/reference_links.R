# List of of reference atlaases links

rl <- data.frame(collection.CSI = c("human",
                                    "human",
                                    "human",
                                    "human",
                                    "mouse",
                                    "mouse",
                                    "mouse"
                                    ),
                 reference.atlas = c("CD4",
                                     "CD8",
                                     "DC",
                                     "MoMac",
                                     "Virus_CD4T",
                                     "Virus_CD8T",
                                     "TILs"
                                     ),
                 name = c("sketched_CD4T_human_ref_v2.rds",
                          "sketched_CD8T_human_ref_v1.rds",
                          "sketched_DC_human_ref_v2.rds",
                          "sketched_MoMac_human_v1.rds",
                          "ref_LCMV_CD4_mouse_release_v1.rds",
                          "ref_CD8_LCMV_mouse_v2.rds",
                          "ref_TILAtlas_mouse_v1.rds"
                          ),
                 figshare_id = c("26310994",
                                 "26310994",
                                 "26310994",
                                 "26310994",
                                 "16592693",
                                 "23764572",
                                 "12478571"
                                 )
)

inst.dir <- "inst/extdata"
dir.create(inst.dir,
           recursive = T)

utils::write.table(rl,
                   file.path(inst.dir, "reference_links.csv"),
                   sep = ",",
                   quote = F,
                   row.names = F)
