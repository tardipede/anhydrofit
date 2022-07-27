library(writexl)
library(anhydrofit)
library(ggplot2)
library(patchwork)

a_Mil = analyse_model(mod_fit_Mil)
a_Par = analyse_model(mod_fit_Par)
a_Ram = analyse_model(mod_fit_Ram)

write_xlsx(a_Mil$pairwise.table, "Milnesium.xlsx")
write_xlsx(a_Par$pairwise.table, "Paramacrobiotus.xlsx")
write_xlsx(a_Ram$pairwise.table, "Ramazzottius.xlsx")



## Mil
p1_Mil = plot_model(mod_fit_Mil)
p2_Mil = plot_ahi(mod_fit_Mil)+theme(axis.text.x = element_text(angle = 90))
p1_Mil/p2_Mil
ggsave(filename = "Milnesium.svg", height = 10, width = 5)
ggsave(filename = "Milnesium.tiff", height = 10, width = 5)

plot_fits(mod_fit_Mil)
ggsave(filename = "Milnesium_predicted.svg", height = 7, width = 7)
ggsave(filename = "Milnesium_predicted.tiff", height = 7, width = 7)

## Par
p1_Par = plot_model(mod_fit_Par)
p2_Par = plot_ahi(mod_fit_Par)+theme(axis.text.x = element_text(angle = 90))
p1_Par/p2_Par
ggsave(filename = "Paramacrobiotus.svg", height = 10, width = 5)
ggsave(filename = "Paramacrobiotus.tiff", height = 10, width = 5)

plot_fits(mod_fit_Par)
ggsave(filename = "Paramacrobiotus_predicted.svg", height = 7, width = 9)
ggsave(filename = "Paramacrobiotus_predicted.tiff", height = 7, width = 9)

## Ram
p1_Ram = plot_model(mod_fit_Ram)
p2_Ram = plot_ahi(mod_fit_Ram)+theme(axis.text.x = element_text(angle = 90))
p1_Ram/p2_Ram
ggsave(filename = "Ramazzottius.svg", height = 10, width = 5)
ggsave(filename = "Ramazzottius.tiff", height = 10, width = 5)

plot_fits(mod_fit_Ram)
ggsave(filename = "Ramazzottius_predicted.svg", height = 5, width = 7)
ggsave(filename = "Ramazzottius_predicted.tiff", height = 5, width = 7)
