rm(list = ls())

library(viridis)
library(ggplot2)

n_eval = 3

#### Individuals ####

n_possibilities = c(50,100,250,500, 1000, 5000, 10000, 12500, 15000, 20000, 25000, 30000, 40000, 50000, 100000, 500000, 10^6)

load("Dossier_AD/Ind_Etape_9_p_100.RData")
list_AD = lists_res_n

load("Dossier_DC/Ind_Etape_14_p_100.RData")
list_DC = lists_res_n

load("Dossier_MT/Ind_Etape_17_p_100.RData")
list_MT = lists_res_n

names(list_AD) = n_possibilities[1:length(list_AD)]
names(list_DC) = n_possibilities[1:length(list_DC)]
names(list_MT) = n_possibilities[1:length(list_MT)]

AD_df = do.call("rbind", lapply(list_AD, data.frame))
AD_df$N = rep(names(list_AD), each = nrow(list_AD[[1]]))
AD_df$time = AD_df$time/10^(9)
AD_df$Method = AD_df$expr

DC_df = do.call("rbind", lapply(list_DC, data.frame))
DC_df$N = rep(names(list_DC), each = nrow(list_DC[[1]]))
DC_df$time = DC_df$time/10^(9)
DC_df$Method = DC_df$expr

MT_df = do.call("rbind", lapply(list_MT, data.frame))
MT_df$N = rep(names(list_MT), each = nrow(list_MT[[1]]))
MT_df$time = MT_df$time/10^(9)
MT_df$Method = MT_df$expr

ggplot_table = rbind(AD_df, DC_df, MT_df)

ggplot_n_ind = aggregate(ggplot_table$time, by =  list(ggplot_table$expr, ggplot_table$N), FUN = function(vect) return(c(mean(as.numeric(vect)), sd(as.numeric(vect)))))
ggplot_n_ind = cbind(ggplot_n_ind[,1:2], ggplot_n_ind$x)
colnames(ggplot_n_ind) = c("Method", "N", "Mean", "SD")

ggplot_n_ind$Compare = rep("Individuals", nrow(ggplot_n_ind))
ggplot_n_ind$Spectral = "Non spectral"
ggplot_n_ind$Spectral[grep("S", ggplot_n_ind$Method)] = "Spectral"
ggplot_n_ind$N = as.numeric(ggplot_n_ind$N)


n_breaks = c(250, 500, 10^3, 5000, 10000, 100000, 500000, 10^6)
n_breaks_label = c(250, 500, 10^3, 5000, "1e+04", 10^5, "5e+05", 10^6)

time_breaks = c(1,2,10,25,50,100,250, 500, 1000,2000)

#### Features ####

p_possibilities = p_possibilities = c(1000, 2500, 5000, 10000, 25000, 50000, 100000, 500000, 10^6)

load("Dossier_AD/Features_Etape_7.RData")
list_AD = lists_res_p

load("Dossier_DC/Features_Etape_7.RData")
list_DC = lists_res_p

load("Dossier_MT/Features_Etape_7.RData")
list_MT = lists_res_p

names(list_AD) = p_possibilities[1:length(list_AD)]
names(list_DC) = p_possibilities[1:length(list_DC)]
names(list_MT) = p_possibilities[1:length(list_MT)]

AD_df = do.call("rbind", lapply(list_AD, data.frame))
AD_df$N = rep(names(list_AD), each = nrow(list_AD[[1]]))
AD_df$time = AD_df$time/10^(9)
AD_df$Method = AD_df$expr

DC_df = do.call("rbind", lapply(list_DC, data.frame))
DC_df$N = rep(names(list_DC), each = nrow(list_DC[[1]]))
DC_df$time = DC_df$time/10^(9)
DC_df$Method = DC_df$expr

MT_df = do.call("rbind", lapply(list_MT, data.frame))
MT_df$N = rep(names(list_MT), each = nrow(list_MT[[1]]))
MT_df$time = MT_df$time/10^(9)
MT_df$Method = MT_df$expr

ggplot_table = rbind(AD_df, DC_df, MT_df)

ggplot_p_feat = aggregate(ggplot_table$time, by =  list(ggplot_table$expr, ggplot_table$N), FUN = function(vect) return(c(mean(as.numeric(vect)), sd(as.numeric(vect)))))
ggplot_p_feat = cbind(ggplot_p_feat[,1:2], ggplot_p_feat$x)
colnames(ggplot_p_feat) = c("Method", "N", "Mean", "SD")

ggplot_p_feat$Compare = rep("Features", nrow(ggplot_p_feat))
ggplot_p_feat$Spectral = "Non spectral"
ggplot_p_feat$Spectral[grep("S", ggplot_p_feat$Method)] = "Spectral"
ggplot_p_feat$N = as.numeric(ggplot_p_feat$N)


# ---- Trying to get the result I want : ----

# Fonctions trouvees : 
# https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/

scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


facet_wrap_custom(~facet_name, scales = "free", ncol = 4, scale_overrides = list(
  scale_override(1, scale_x_continuous(breaks = c(5750, 5900))),
  scale_override(6, scale_x_continuous(breaks = c(17800, 17900)))
))


# ---- Figure : ----

ggplot_all = rbind(ggplot_n_ind, ggplot_p_feat)

ggplot_all$Spec_comp = paste(ggplot_all$Compare, ggplot_all$Spectral, sep = " - ")

# Changement noms de methodes pour coller avec le reste de l'article
Method_vect = as.character(ggplot_all$Method)
Method_vect[which(ggplot_all$Method == "ACuni")] = "AD" 
Method_vect[which(ggplot_all$Method == "MTuni")] = "MC" 
Method_vect[which(ggplot_all$Method == "SDC")] = "SpDC" 
Method_vect[which(ggplot_all$Method == "SACuni")] = "SpAD" 
Method_vect[which(ggplot_all$Method == "SMTuni")] = "SpMC" 

ggplot_all$Method = factor(Method_vect)
ggplot_all = ggplot_all[which(ggplot_all$N>=250),]

breaks_var = c(1000, 2500, 5000, 10000, 25000, 50000, 10^5, 10^6)
labels_var =  c(1000, 2500, 5000, 10000, 25000, 50000, "1e+05", "1e+06")

breaks_ind = c(250, 500, 2500, 10^4,  10^5, 10^6)
labels_ind = c(250, 500, 2500, "1e+04","1e+05", "1e+06")

# Trace de la figure
gg = ggplot(data = ggplot_all, aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  theme_bw() + ylab("Time (s)") + xlab("") + 
  theme(legend.position="bottom", legend.box = "horizontal") +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_viridis(discrete = TRUE) +
  facet_wrap_custom( ~ Spec_comp, ncol = 2, scales = "free", scale_overrides = list(
    scale_override(1, scale_x_continuous(trans = "log", breaks = breaks_var, labels = labels_var)), 
    scale_override(1, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500, 2500))),
    
    scale_override(2, scale_x_continuous(trans = "log", breaks = breaks_var, labels = labels_var)), 
    scale_override(2, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500, 2500))),
    
    scale_override(3, scale_x_continuous(trans = "log", breaks = breaks_ind, labels = labels_ind)),
    scale_override(3, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500, 2500))),
    
    scale_override(4, scale_x_continuous(trans = "log", breaks = breaks_ind, labels = labels_ind)),
    scale_override(4, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500)))
    
  ))
gg

ggsave(plot = gg, filename = "Figure2.pdf", device = "pdf", width = 25, height = 15, units = "cm")
