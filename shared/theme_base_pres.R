library(ggplot2)
theme_set(theme_minimal(base_size = 18)+
            theme(text = element_text(size=12),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  axis.title.x = element_text(size=16,margin = margin(t = -2)),
                  axis.title.y = element_text(size=16),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom',
                  legend.box.spacing = unit(0,'pt'),
                  legend.key.spacing.y = unit(-8, 'pt'),
                  legend.margin = margin(t=0,r=0,b=0,l=0),
                  legend.box.margin = margin(t=0,r=0,b=0,l=0),
                  legend.text = element_text(size=18,margin=margin(0,0,0,0)),
                  plot.margin = margin(t=4,l=4,r=0,b=0))
)
