library(ggplot2)

# データの準備 ----------------------------------------------------

# 2. 関心のあるパスウェイを指定
selected_pathways <- c(
  "Warburg Effect",
  "Gluconeogenesis",
  "Glycolysis",
  "Citric Acid Cycle",
  "Pyruvaldehyde Degradation",
  "Urea Cycle",
  "Glutamate Metabolism",
  "Glycine and Serine Metabolism",
  "Pyruvate Metabolism",
  "Spermidine and Spermine Biosynthesis"
)

# 全範囲
df0 <- as.data.frame(B2$`Range of p-values`)
df0 <- df0[row.names(df0) %in% selected_pathways, ]
df0$pathway <- row.names(df0)
colnames(df0)[1:3] <- c("p_lower", "p_median", "p_upper")
df0$method <- "All range"  # ラベル付け

# B4: 通常の方法
#df1 <- as.data.frame(B4$`Range of p-values for each pathway (95% confidence interval)`)
df1 <- as.data.frame(B4[[1]])
df1 <- df1[row.names(df1) %in% selected_pathways, ]
df1$pathway <- row.names(df1)
colnames(df1)[1:3] <- c("p_lower", "p_median", "p_upper")
df1$method <- "Binomial CI (95%)"  # ラベル付け

# B5: Beta分布の方法
#df2 <- as.data.frame(B5$`Range of p-values for each pathway (95% CI using Beta distribution)`)
df2 <- as.data.frame(B5[[1]])
df2 <- df2[row.names(df2) %in% selected_pathways, ]
df2$pathway <- row.names(df2)
colnames(df2)[1:3] <- c("p_lower", "p_median", "p_upper")
df2$method <- "Binomial CI (80%)"  # ラベル付け

# 統合して因子順序を揃える ---------------------------------------

df_all <- rbind(df0, df1, df2)
df_all$pathway <- factor(df_all$pathway, levels = rev(selected_pathways))
df_all$method <- factor(df_all$method, levels = c("All range", "Binomial CI (95%)", "Binomial CI (80%)"))

# プロット -------------------------------------------------------

# 軸の範囲を明示的に設定
y_min <- min(df_all$p_lower)
y_max <- max(df_all$p_upper)

p <- ggplot(df_all, aes(x = pathway, y = p_median)) +
  geom_point(color = "blue") +
  geom_errorbar(aes(ymin = p_lower, ymax = p_upper), width = 0.2) +
  coord_flip() +
  scale_y_log10(
    limits = c(y_min, y_max),
    breaks = c(0.05)
  ) +
  facet_wrap(~ method, ncol = 3) +
  ylab("p-value") +
  xlab("Pathway") +
  ggtitle("Comparison of p-value 95% CIs by Method") +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  )

print(p)

