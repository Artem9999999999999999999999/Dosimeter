# Загрузка необходимых библиотек
library(DESeq2)
library(ggplot2)
library(dplyr)

# Функция для загрузки аннотации и фильтрации данных
load_and_filter_annotation <- function(annotation_file) {
  # Загрузка аннотации
  annotation <- read.gff(annotation_file, na.strings = c(".", "?"), GFF3 = FALSE)
  
  # Фильтрация аннотации по генам
  annotation_filter <- annotation %>%
    filter(feature == "gene") %>%
    dplyr::select(c(1, 3, 7, 9)) %>% 
    setNames(c("chr", "feature", "strand", "attr")) %>% 
    mutate(id = sub('.*gene_id "([^"]+)".*', '\\1', attr),
           gene_name = sub('.*gene_symbol "([^"]+)".*', '\\1', attr),
           chr = paste0("chr", chr)) %>% 
    dplyr::select(-attr)
  
  return(annotation_filter)
}

# Загрузка и фильтрация аннотации
annotation_file <- "data_r6.55_bw2_RS/dmel-all-r6.55.gtf.gz"
annotation_filter <- load_and_filter_annotation(annotation_file)

# Функция для объединения результатов DE анализа с аннотацией
merge_annotation_and_results <- function(annotation, diff_expr_results, file_path) {
  # Объединение результатов DE анализа с аннотацией
  merged_results <- merge(select(annotation, id, chr, strand, gene_name), diff_expr_results, by = "id", all.y = TRUE)
  
  # Сохранение объединенных результатов в файл
  write.csv(merged_results, file = file_path, row.names = TRUE)
}

# Объединение результатов DE анализа с аннотацией для CHD1
merge_annotation_and_results(annotation_filter, diff_expr_CHD1, "results_r6.55_bw2_RS/diff_expr_CHD1(r6.55)bw2_RS.csv")

# Фильтрация результатов по хромосомам и создание группировки
filtered_counts_CHD1.df <- subset(merged_results, !(chr %in% c("chrY", "chrmitochondrion_genome", "chrrDNA")))
filtered_counts_CHD1.df$group <- ifelse(filtered_counts_CHD1.df$chr == "chrX", "chrX", "Autosome")

# Разделение результатов на повышенную и пониженную экспрессию
filtered_counts_CHD1_up.df <- filtered_counts_CHD1.df[filtered_counts_CHD1.df$log2FoldChange > 0, ]
filtered_counts_CHD1_down.df <- filtered_counts_CHD1.df[filtered_counts_CHD1.df$log2FoldChange < 0, ]

# Создание папки для результатов
results_folder <- "results_r6.55_bw2_RS"
if (!file.exists(results_folder)) {
  dir.create(results_folder)
}

# Сохранение результатов в файлы
write.csv(filtered_counts_CHD1_up.df, file.path(results_folder, "CHD1_UP(r6.55)bw2.RS.csv"), row.names = FALSE)
write.csv(filtered_counts_CHD1_down.df, file.path(results_folder, "CHD1_DOWN(r6.55)bw2.RS.csv"), row.names = FALSE)

# Создание графиков
color_palette <- c("#92c5de", "#f4a582")

# Создание boxplot для групп с повышенной экспрессией
ggplot(filtered_counts_CHD1_up.df, aes(x = group, y = log2FoldChange, fill = group)) +
  geom_boxplot() +
  labs(title = "Boxplot of log2FC for DE genes of PRO-Seq UP (r6.55)", x = "Group", y = "log2FoldChange") +
  scale_fill_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  ylim(0, 1)

# Создание boxplot для групп с пониженной экспрессией
ggplot(filtered_counts_CHD1_down.df, aes(x = group, y = log2FoldChange, fill = group)) +
  geom_boxplot() +
  labs(title = "Boxplot of log2FC for DE genes of PRO-Seq DOWN (r6.55)", x = "Group", y = "log2FoldChange") +
  scale_fill_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  ylim(-1, 0)

# Создание density plot для групп с повышенной экспрессией
ggplot(filtered_counts_CHD1_up.df, aes(x = log2FoldChange, color = group)) +
  geom_density(alpha = 0, size = 1) + 
  labs(title = "Density Plot of log2FC for DE genes of PRO-Seq UP (r6.55)", x = "log2FoldChange", y = "Density") +
  scale_color_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  xlim(0, 2)

# Создание density plot для групп с пониженной экспрессией
ggplot(filtered_counts_CHD1_down.df, aes(x = log2FoldChange, color = group)) +
  geom_density(alpha = 0, size = 1) + 
  labs(title = "Density Plot of log2FC for DE genes of PRO-Seq DOWN (r6.55)", x = "log2FoldChange", y = "Density") +
  scale_color_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  xlim(-2, 0)

# Выполнение статистических тестов
t_test_result_up <- t.test(log2FoldChange ~ group, data = filtered_counts_CHD1_up.df)
t_test_result_down <- t.test(log2FoldChange ~ group, data = filtered_counts_CHD1_down.df)

# Вывод результатов статистических тестов
print("T-Test UP:")
print(t_test_result_up)
print("T-Test DOWN:")
print(t_test_result_down)

# Выполнение рангового теста Манна-Уитни
wilcox_test_result_up <- wilcox.test(log2FoldChange ~ group, data = filtered_counts_CHD1_up.df)
wilcox_test_result_down <- wilcox.test(log2FoldChange ~ group, data = filtered_counts_CHD1_down.df)

# Вывод результатов рангового теста Манна-Уитни
print("Mann-Whitney U Test UP:")
print(wilcox_test_result_up)

print("Mann-Whitney U Test DOWN:")
print(wilcox_test_result_down)

# Вычисление медиан
medians_up <- with(filtered_counts_CHD1_up.df, tapply(log2FoldChange, group, median))
medians_down <- with(filtered_counts_CHD1_down.df, tapply(log2FoldChange, group, median))

# Вывод медиан
cat("\nMedians UP:\n")
print(medians_up)
cat("\nMedians DOWN:\n")
print(medians_down)

# Создание папки для сохранения графиков
figures_folder <- "figures_r6.55_bw2_RS"
if (!file.exists(figures_folder)) {
  dir.create(figures_folder)
}

# Сохранение графиков в файлы
boxplot_path_up <- file.path(figures_folder, "Boxplot_DE_chd1_up(r6.55)bw2RS.png")
boxplot_path_down <- file.path(figures_folder, "Boxplot_DE_chd1_down(r6.55)bw2RS.png")
density_plot_path_up <- file.path(figures_folder, "Density_plot_DE_chd1_up(r6.55)bw2RS.png")
density_plot_path_down <- file.path(figures_folder, "Density_plot_DE_chd1_down(r6.55)bw2RS.png")

png(boxplot_path_up)
ggplot(filtered_counts_CHD1_up.df, aes(x = group, y = log2FoldChange, fill = group)) +
  geom_boxplot() +
  labs(title = "Boxplot of log2FC for DE genes of Chd1 UP(r6.55) bw2 RS", x = "Group", y = "log2FoldChange") +
  scale_fill_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  ylim(0, 1)
dev.off()

png(boxplot_path_down)
ggplot(filtered_counts_CHD1_down.df, aes(x = group, y = log2FoldChange, fill = group)) +
  geom_boxplot() +
  labs(title = "Boxplot of log2FC for DE genes of Chd1 DOWN(r6.55) bw2 RS", x = "Group", y = "log2FoldChange") +
  scale_fill_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  ylim(-1, 0)
dev.off()

png(density_plot_path_up)
ggplot(filtered_counts_CHD1_up.df, aes(x = log2FoldChange, color = group)) +
  geom_density(alpha = 0, size = 1) + 
  labs(title = "Density Plot of log2FC for DE genes of Chd1 UP(r6.55) bw2 RS", x = "log2FoldChange", y = "Density") +
  scale_color_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  xlim(0, 2)
dev.off()

png(density_plot_path_down)
ggplot(filtered_counts_CHD1_down.df, aes(x = log2FoldChange, color = group)) +
  geom_density(alpha = 0, size = 1) + 
  labs(title = "Density Plot of log2FC for DE genes of Chd1 DOWN(r6.55) bw2 RS", x = "log2FoldChange", y = "Density") +
  scale_color_manual(values = color_palette, name = "Group") +
  theme_minimal() +
  xlim(-2, 0)
dev.off()

cat("Графики сохранены в папке:", figures_folder, "\n")
