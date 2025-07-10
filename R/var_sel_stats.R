# multi approach variable selection script

source('xrf_functions.R')

library(tidyverse)
library(corrplot)
library(cluster)
library(factoextra)
library(randomForest)
library(Boruta)
library(fpc)



# combined contains all the variables in one data frame
combined <- combined_df_creator()

# variable selection
sel_vars <- c("Fe", "Mn", "S", "Ti", "Al", "K", "Si", "Ca", "Ba", "Zr", 
              "Mn_Fe", "Fe_Ti", "S_Ti", 'Si_Ti', "Ba_Ti", "Ca_Ti", "Ti_Al", "Zr_Rb", 
              "TCWt", 'TNWt', 'TSWt', "RABD673", "Moving.Average")

df_sel <- combined %>%
  select(any_of(sel_vars)) %>%
  drop_na()

#### Hierarchical clustering ##################################################
# corr matrix plot
cor_mat <- cor(df_sel, use = "pairwise.complete.obs")
corrplot(cor_mat, method = "color", type = "upper", tl.cex = 0.7)

combined %>% 
  select(any_of(sel_vars)) %>%
  create_corr_map()


hc <- hclust(as.dist(1 - abs(cor_mat)), method = "ward.D2")
plot(hc, main = "Variable Clustering Dendrogram")

# Group into clusters (adjust k)
groups <- cutree(hc, k = 5)
group_df <- data.frame(variable = names(groups), group = groups)
print(group_df %>% arrange(group))



#### Random forest ############################################################

rf <- randomForest(df_sel, proximity = TRUE, importance = TRUE)
varImpPlot(rf, main = "Unsupervised Random Forest Variable Importance")
importance_df <- importance(rf) %>% as.data.frame() %>%
  rownames_to_column("variable") %>%
  arrange(desc(MeanDecreaseGini))
print(importance_df)

#### Boruta  ##################################################################
set.seed(123)
df_sel_wPos <- df_sel %>%
  mutate(position..mm. = facies_sync$position..mm.)
df_boruta <- df_sel_wPos %>%
  left_join(facies_sync, by = "position..mm.")

boruta_out <- Boruta(x = select(df_boruta, -c(facies_class, position..mm.)),
                     y = as.factor(df_boruta$facies_class),
                     doTrace = 1, maxRuns = 100)

plot(boruta_out, las = 2, cex.axis = 0.7)
getSelectedAttributes(boruta_out, withTentative = TRUE)
#### Cluster stability with variable sets #####################################
set.seed(123)

cluststab <- clusterboot(
  scale(df_sel),
  B = 100,
  clustermethod = kmeansCBI,
  krange = 2:5,
  seed = 123,
  count = TRUE
)

print(cluststab$bootmean)
print(cluststab$result$partition)


#### PCA var explained ########################################################
pca_out <- prcomp(df_sel, scale. = TRUE)
fviz_eig(pca_out, addlabels = TRUE, ylim = c(0, 50))

summary(pca_out)$importance[2:3, 1:10]  # First 10 components
