rm(list = ls())

# library(viridis)
library(ggplot2)

n_eval = 3

# setwd("/travail/ahulot/Simulations_timings_2020_01_15/Dossier_res/")
setwd("/home/hulot/Documents/packages_R/mergeTrees/inst/article_reviews_2020_01_10/Simulations_timing_2020_01_15/Dossier_res/")

load("Features_Etape_9_n_100_q_3.RData")
load("Ind_Etape_13_p_100_q_3.RData")
load("Tables_Etape_13_n_100_p_100.RData")

names(lists_res_n)
names(lists_res_p)
names(lists_res_q)

# p_possibilities = c(1000, 2500, 5000, 10000, 25000, 50000, 100000, 500000, 10^6)
# names(lists_res_p)=p_possibilities[1:length(lists_res_p)]


Feat_df = do.call("rbind", lapply(lists_res_p, data.frame))
Feat_df$N = rep(names(lists_res_p), each = nrow(lists_res_p[[1]]))
Feat_df$time = Feat_df$time/10^(9)
Feat_df$Comp = "Features"

Ind_df = do.call("rbind", lapply(lists_res_n, data.frame))
Ind_df$N = rep(names(lists_res_n), each = nrow(lists_res_n[[1]]))
Ind_df$time = Ind_df$time/10^(9)
Ind_df$Comp = "Individuals"

Tab_df = do.call("rbind", lapply(lists_res_q, data.frame))
Tab_df$N = rep(names(lists_res_q), each = nrow(lists_res_q[[1]]))
Tab_df$time = Tab_df$time/10^(9)
Tab_df$Comp = "Tables"

ggplot_table = rbind(Tab_df, Ind_df, Feat_df)

ggplot_mean = aggregate(ggplot_table$time, by =  list(ggplot_table$expr, ggplot_table$N, ggplot_table$Comp), 
                        FUN = function(vect) return(c(mean(as.numeric(vect)), sd(as.numeric(vect)))))
ggplot_mean = cbind(ggplot_mean[,1:3], ggplot_mean$x)
colnames(ggplot_mean) = c("Method", "N", "Comparison", "Mean", "SD")
ggplot_mean$N = as.numeric(as.character(ggplot_mean$N))
ggplot_mean$Spectral = "No Spectral"
ggplot_mean$Spectral[grep("sp", ggplot_mean$Method)] = "Spectral"

gg = ggplot(data = ggplot_mean, aes(x = N, y = Mean, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  scale_x_continuous(trans = "log", breaks = c(100, 500, 1000, 10^4, 10^5, 10^6)) +
  scale_y_continuous(trans = "log", breaks = c(0.01,0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 100, 10^3, 2500, 5000)) +
  theme_bw() + ylab("Time (s)") + xlab("") + 
  theme(legend.position="bottom", legend.box = "horizontal", strip.placement = "outside",
        strip.background.x = element_blank()) +
  guides(color = guide_legend(nrow = 1)) +
  facet_grid(Spectral~Comparison, scales = "free", switch = "x")
gg

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
  # facet_super <- facet_grid(...)
  
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


# facet_wrap_custom(~facet_name, scales = "free", ncol = 4, scale_overrides = list(
# scale_override(1, scale_x_continuous(breaks = c(5750, 5900))),
# scale_override(6, scale_x_continuous(breaks = c(17800, 17900)))
# ))


# ---- Figure : ----

breaks_var = c(1000, 2500, 5000, 10000, 25000, 50000, 10^5, 10^6)
labels_var =  c(1000, 2500, 5000, 10000, 25000, 50000, "1e+05", "1e+06")

breaks_ind = c(250, 500, 2500, 10^4,  10^5, 10^6)
labels_ind = c(250, 500, 2500, "1e+04","1e+05", "1e+06")

breaks_tab = c(250, 500, 2500, 10^4,  10^5, 10^6)
labels_tab = c(250, 500, 2500, "1e+04","1e+05", "1e+06")

ggplot_m

gg = ggplot(data = ggplot_mean, aes(x = N, y = Mean, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  scale_x_continuous(trans = "log", breaks = c(100, 500, 1000, 10^4, 10^5, 10^6)) +
  scale_y_continuous(trans = "log", breaks = c(0.01,0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 100, 10^3, 2500, 5000)) +
  theme_bw() + ylab("Time (s)") + xlab("") + 
  theme(legend.position="bottom", legend.box = "horizontal", strip.placement = "outside",
        strip.background.x = element_blank()) +
  guides(color = guide_legend(nrow = 1)) +
  facet_wrap_custom( Spectral~Comparison,
                     ncol = 3, scales = "free",
                     strip.position = c("bottom"),
                     scale_overrides = list(
                       scale_override(1, scale_x_continuous(trans = "log", breaks = breaks_ind, labels = labels_ind)),
                       scale_override(1, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500, 2500))),
                       
                       scale_override(2, scale_x_continuous(trans = "log", breaks = breaks_ind, labels = labels_ind)),
                       scale_override(2, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500)))
                       
                     )) 
gg

gg1 = ggplot(data = ggplot_mean[ggplot_mean$Comparison=="Individuals",], aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  theme_bw() + ylab("Time (s)") +
  theme(legend.position="bottom", legend.box = "horizontal", strip.background = element_blank()) +
  guides(color = guide_legend(nrow = 1)) +
  # scale_color_viridis(discrete = TRUE) +
  facet_wrap_custom(Spectral~.,
                    # Compare ~ Spectral,
                    ncol = 1, scales = "free",
                    strip.position = c("top"),
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(trans = "log", breaks = breaks_ind, labels = labels_ind)),
                      scale_override(1, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500, 2500))),

                      scale_override(2, scale_x_continuous(trans = "log", breaks = breaks_ind, labels = labels_ind)),
                      scale_override(2, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500)))

                    )) +
  xlab("Number of Individuals")
gg1

gg2 = ggplot(data = ggplot_mean[ggplot_mean$Comparison=="Features",], aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  theme_bw() + ylab("Time (s)") +
  theme(legend.position="bottom", legend.box = "horizontal", strip.background = element_blank()) +
  guides(color = guide_legend(nrow = 1)) +
  # scale_color_viridis(discrete = TRUE) +
  facet_wrap_custom(~Spectral,
                    # Compare ~ Spectral,
                    ncol = 1, scales = "free",
                    strip.position = c("top"),
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(trans = "log", breaks = breaks_var, labels = labels_var)),
                      scale_override(1, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500, 2500))),

                      scale_override(2, scale_x_continuous(trans = "log", breaks = breaks_var, labels = labels_var)),
                      scale_override(2, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500)))

                    )) +
  xlab("Number of Features/Trees")
gg2

gg3 = ggplot(data = ggplot_mean[ggplot_mean$Comparison=="Tables",], aes(x = N, y = Mean, colour = Method)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=Mean-sqrt(1.96)*SD/sqrt(n_eval), ymax=Mean+sqrt(1.96)*SD/sqrt(n_eval)), width=.1) +
  scale_x_continuous(trans = "log")+
  scale_y_continuous(trans = "log")+
  theme_bw() + ylab("Time (s)") +
  guides(color = guide_legend(nrow = 1)) +
  # scale_color_viridis(discrete = TRUE) +
  facet_wrap_custom(~Spectral, 
                    ncol = 1, scales = "free",
                    strip.position = c("top"),
                    scale_overrides = list(
                      scale_override(1, scale_x_continuous(trans = "log", breaks = breaks_tab, labels = labels_tab)),
                      scale_override(1, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500, 2500))),
                      
                      scale_override(2, scale_x_continuous(trans = "log", breaks = breaks_tab, labels = labels_tab)),
                      scale_override(2, scale_y_continuous(trans = "log", breaks = c(0.05, 0.25, 1, 5, 25,  100, 500)))
                      
                    )) +
  theme(legend.position="bottom", legend.box = "horizontal", strip.background = element_blank()) +
  xlab("Number of datasets")
gg3

library(ggpubr)
gg = ggarrange(gg1, gg2, gg3,  common.legend = TRUE, legend="bottom", ncol = 3)
gg
ggsave(gg, filename = "Figure2.pdf", device = "pdf", width = 25, height = 15, units = "cm")

