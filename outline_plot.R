library(ggplot2)
library(patchwork)
library(dplyr)
library(tiff)

selected <- data[grepl(
  "5351485e-a9b7-4d61-9f9d-f7a6737dd6ec",
  data$CellID
), ]

########################################
# LANDMARKS
########################################
com <- as.numeric(strsplit(selected$Centre_of_mass, " - ")[[1]])

cx <- com[1]
cy <- com[2]


img_h <- dim(img)[1]
img_w <- dim(img)[2]
landmarks <- data.frame(
  x = as.numeric(selected[, grepl("RawCoordinates_X", names(selected))]),
  y = as.numeric(selected[, grepl("RawCoordinates_Y", names(selected))])
) %>%
  filter(!is.na(x), !is.na(y))

landmarks$Position <- seq_len(nrow(landmarks))
landmarks$x_crop <- landmarks$x - cx
landmarks$y_crop <- landmarks$y - cy

x_rot <- -landmarks$y_crop
y_rot <- landmarks$x_crop

landmarks$x_plot <- x_rot + img_w/2
landmarks$y_plot <- y_rot + img_h/2

landmarks_closed <- rbind(
  landmarks,
  landmarks[1, ]
)

landmarks_closed$Position <- c(
  landmarks$Position,
  landmarks$Position[1]
)
outline_pts <- landmarks %>%
  filter(Position %% 5 == 0)

p_outline <- ggplot(landmarks_closed,
                    aes(x_plot,y_plot)) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 0.8) +
  coord_equal() +
  theme_bw() +
  ggtitle("Outline Landmarks")
p_outline
########################################
# IMAGE
########################################

img <- readTIFF("D:/2024-2025_sperm_images/Deer_Mouse/Combined_relabelled_deer_mice/2025-01-27_13-21-36/Normal/Slide_20_150X_L10_tile_x005_y007.tif-57875_5351485e-a9b7-4d61-9f9d-f7a6737dd6ec.tiff")

img_df <- expand.grid(
  x = seq_len(dim(img)[2]),
  y = seq_len(dim(img)[1])
)

img_df$intensity <- as.vector(img[nrow(img):1, ])

library(ggnewscale)

p_image <- ggplot() +

  geom_raster(
    data = img_df,
    aes(x, y, fill = intensity)
  ) +

  scale_fill_gradient(
    low = "black",
    high = "white"
  ) +

  new_scale_colour() +

  geom_path(
    data = landmarks_closed,
    aes(
      x_plot,
      y_plot,
      colour = Position
    ),
    linewidth = 2
  ) +

  geom_point(
    data = landmarks_closed,
    aes(
      x_plot,
      y_plot,
      colour = Position
    ),
    size = 0.5
  ) +
  geom_point(
    data = outline_pts,
    aes(
      x_plot,
      y_plot,
      colour = Position
    ),
    size = 2
  )+

  scale_colour_viridis_c(
    option = "turbo",
    limits = c(1,100)
  ) +

  coord_equal() +
  theme_void() +
  ggtitle("Cell and Outline")+
  guides(
    fill = "none"
  )
p_image

########################################
# PROFILES
########################################

angle_profile <- as.numeric(
  data[i, grep("^Angle", names(data))]
)

radius_profile <- as.numeric(
  data[i, grep("^Radius", names(data))]
)

diameter_profile <- as.numeric(
  data[i, grep("^Diameter", names(data))]
)

angle_df <- data.frame(
  Position = seq_along(angle_profile),
  Value = angle_profile
)

radius_df <- data.frame(
  Position = seq_along(radius_profile),
  Value = radius_profile
)

diameter_df <- data.frame(
  Position = seq_along(diameter_profile),
  Value = diameter_profile
)
#######

#######
profile_plot <- function(df, ylab){

  seg_df <- data.frame(
    x = df$Position[-nrow(df)],
    xend = df$Position[-1],
    y = df$Value[-nrow(df)],
    yend = df$Value[-1],
    Position = df$Position[-nrow(df)]
  )

  ggplot(seg_df) +
    geom_segment(
      aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend,
        colour = Position
      ),
      linewidth = 1
    ) +
    scale_colour_viridis_c(
      option = "turbo",
      limits = c(1, 100)
    ) +
    theme_bw() +
    xlab("Profile Position") +
    ylab(ylab)+
    theme(
      legend.position = "none"
    )
}

########################################
# ANGLE GRAPH
########################################

p_angle <- profile_plot(df = angle_df,ylab = "Angle")

########################################
# RADIUS GRAPH
########################################

p_radius <- profile_plot(radius_df, ylab = "Radius")

########################################
# DIAMETER GRAPH
########################################

p_diameter <- profile_plot(diameter_df, ylab = "Diameter")
########################################
# PATCHWORK
########################################

final_plot <-
  (p_image |p_angle  / p_radius / p_diameter)

final_plot

ggsave(
  "Cell_Profile_Panel.pdf",
  final_plot,
  width = 12,
  height = 8
)
