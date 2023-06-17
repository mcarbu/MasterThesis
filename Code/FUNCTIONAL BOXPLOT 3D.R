## FUNCTIONAL BOXPLOT 2.5D ## 

# Libraries
library(fdaPDE)
library(plotly)
library(crosstalk)
library(roahd)

options(warn = -1)

# Useful functions
volume_plot <- function(FEM, ...){
  plot_ly(type = 'mesh3d',
          
          x = FEM$FEMbasis$mesh$nodes[,1],
          y = FEM$FEMbasis$mesh$nodes[,2],
          z = FEM$FEMbasis$mesh$nodes[,3],
          
          i = FEM$FEMbasis$mesh$faces[,1]-1,
          j = FEM$FEMbasis$mesh$faces[,2]-1,
          k = FEM$FEMbasis$mesh$faces[,3]-1, 
          flatshading = T,
          
          intensity = FEM$coeff,
          color = FEM$coeff,
          ...) %>%
    colorbar(len = 1)
}
plot_functional_boxplot <- function(FEMbasis, lowerfence, q1, median, q3, upperfence, colorscale = 'Set3', ...){
  
  ## Layout options 
  axx <- list(
    title = '',
    showgrid = F,
    zeroline = F,
    showticklabels = F 
  )
  
  add_ann <- function(text){
    return(list(text = text,
                xref = "paper",
                yref = "paper",
                yanchor = "bottom",
                xanchor = "center",
                align = "center",
                x = 0.5,
                y = 0.8,
                showarrow = FALSE))
  }
  
  set_camera <- list(eye = list(x = 1.8,y = 1.8, z = 1.8))
  
  ## Univariate boxplot:integral over the domain
  box <- plot_ly(x = 0,
                 type = "box", 
                 lowerfence=mean(lowerfence),
                 q1=mean(q1), 
                 median=mean(median),
                 q3=mean(q3),
                 upperfence=mean(upperfence), 
                 quartilemethod="inclusive",
                 showlegend = F,
                 color = I('darkgrey')) %>%
    layout(xaxis = axx,
           yaxis = axx)
  
  
  ## Summary statistics hubs 
  p_me <- volume_plot(FEM(median, FEMbasis), scene = 'scene3', coloraxis = 'coloraxis') %>%
    layout(annotations = add_ann('Median'))
  p_q1 <- volume_plot(FEM(q1, FEMbasis), scene = 'scene2', coloraxis = 'coloraxis') %>%
    layout(annotations = add_ann('Q1'))
  p_q3 <- volume_plot(FEM(q3, FEMbasis), scene = 'scene4', coloraxis = 'coloraxis') %>%
    layout(annotations = add_ann('Q3'))
  p_lf <- volume_plot(FEM(lowerfence, FEMbasis), coloraxis = 'coloraxis') %>%
    layout(annotations = add_ann('Lower fence'))
  p_uf <- volume_plot(FEM(upperfence, FEMbasis), scene = 'scene5', coloraxis = 'coloraxis') %>%
    layout(annotations = add_ann('Upperfence'))
  
  
  ## Subplots
  fig <- subplot(subplot(p_lf,p_q1,p_me,p_q3,p_uf),box,nrows = 2, titleX = T, titleY = T) %>% 
    layout(title = "Functional boxplot",
           scene = list(domain=list(x=c(0,0.2),y=c(0.5,1)),
                        xaxis=axx, yaxis=axx, zaxis=axx,
                        aspectmode='cube', camera = set_camera),
           scene2 = list(domain=list(x=c(0.2,0.4),y=c(0.5,1)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         aspectmode='cube', camera = set_camera),
           scene3 = list(domain=list(x=c(0.4,0.6),y=c(0.5,1)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         aspectmode='cube', camera = set_camera),
           scene4 = list(domain=list(x=c(0.6,0.8),y=c(0.5,1)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         aspectmode='cube', camera = set_camera), 
           scene5 = list(domain=list(x=c(0.8,1),y=c(0.5,1)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         aspectmode='cube', camera = set_camera),
           coloraxis=list(colorscale=colorscale, 
                          cmin = min(c(lowerfence,q3,median,q1,upperfence)),
                          cmax = max(c(lowerfence,q3,median,q1,upperfence)), 
                          colorbar = list(nticks = 10, 
                                          tickmode = 'auto',
                                          len = 0.4,
                                          y = .75))
    )
  fig
}
plot_functional_boxplot_data <- function(FEMbasis, lowerfence, q1, median, q3, upperfence, data, colorscale = 'Set3', ...){
  # Number of functions
  N <- dim(data)[2]
  
  # Set axis layout
  axx <- list(
    title = '',
    showgrid = F,
    zeroline = F,
    showticklabels = F 
  )
  
  # Set hovertemplate of the functional boxplot
  hovertemp <- paste('x: %{x:.2f}',
                     '<br>y: %{y:.2f}',
                     '<br>z: %{z:.2f}',
                     '<br><b>Value</b>: %{intensity:.2f}')
  
  # Titles of functional boxplot objects
  add_ann <- function(text, x){
    return(list(text = text,
                xref = "paper",
                yref = "paper",
                yanchor = "bottom",
                xanchor = "center",
                align = "center",
                x = x,
                y = 0.9,
                showarrow = FALSE))
  }
  
  # Set Opacity
  set_opacity <- function(fun){
    a <- sum(fun <= lowerfence)
    b <- sum(fun >= lowerfence & fun <= q1)
    c <- sum(fun >= q1 & fun <= median)
    d <- sum(fun >= median & fun <= q3)
    e <- sum(fun >= q3 & fun <= upperfence)
    f <- sum(fun >= upperfence)
    prevalence <- which.max(c(a,b,c,d,e,f))
    if(prevalence == 1)
      return(as.list(c(1,.5,.5,.5,.5,rep(1,N))))
    else if(prevalence == 2)
      return(as.list(c(1,1,.5,.5,.5,rep(1,N))))
    else if(prevalence == 3)
      return(as.list(c(.5,1,1,.5,.5,rep(1,N))))
    else if(prevalence == 4)
      return(as.list(c(.5,.5,1,1,.5,rep(1,N))))
    else if(prevalence == 5)
      return(as.list(c(.5,.5,.5,1,1,rep(1,N))))
    else
      return(as.list(c(.5,.5,.5,.5,1,rep(1,N))))
  }
  
  # Set camera perspective
  set_camera <- list(eye = list(x = 2,y = 2, z = 2))
  
  # Plot Summary statistics
  p_me <- volume_plot(FEM(median, FEMbasis), scene = 'scene3', coloraxis = 'coloraxis',
                       hovertemplate = hovertemp) %>%
    layout(annotations = add_ann('Median', 1.2))
  p_q1 <- volume_plot(FEM(q1, FEMbasis), scene = 'scene2', coloraxis = 'coloraxis',
                       hovertemplate = hovertemp) %>%
    layout(annotations = add_ann('Q1', 0.9))
  p_q3 <- volume_plot(FEM(q3, FEMbasis), scene = 'scene4', coloraxis = 'coloraxis',
                       hovertemplate = hovertemp) %>%
    layout(annotations = add_ann('Q3', 1.5))
  p_lf <- volume_plot(FEM(lf, FEMbasis), coloraxis = 'coloraxis',
                       hovertemplate = hovertemp) %>%
    layout(annotations = add_ann('Lower fence', 0.7))
  p_uf <- volume_plot(FEM(uf, FEMbasis), scene = 'scene5', coloraxis = 'coloraxis',
                       hovertemplate = hovertemp) %>%
    layout(annotations = add_ann('Upper fence', 1.7))
  
  # Default button
  default_v <- c(rep(T,5),rep(F,N))
  default_op <- c(rep(1,5),rep(1,N))
  button = list()
  button[[1]] = list(method = "restyle",
                     args = list("visible", as.list(default_v)),
                     args2 = list("opacity", as.list(default_op)),
                     label = "None")
  
  # Plot NON-visible traces
  p <- plot_ly()
  for (i in 1:N) {
    p <- p %>%
      add_trace(type = 'mesh3d',
                
                x = FEMbasis$mesh$nodes[,1],
                y = FEMbasis$mesh$nodes[,2],
                z = FEMbasis$mesh$nodes[,3],
                
                i = FEMbasis$mesh$faces[,1]-1,
                j = FEMbasis$mesh$faces[,2]-1,
                k = FEMbasis$mesh$faces[,3]-1, 
                flatshading = T,
                
                intensity = data[,i],
                color = data[,i],
                scene = 'scene6', 
                coloraxis = 'coloraxis', 
                visible = F,
                hovertemplate = paste('x: %{x:.2f}',
                                      '<br>y: %{y:.2f}',
                                      '<br>z: %{z:.2f}',
                                      '<br><b>Value</b>: %{intensity:.2f}', 
                                      '<br>Data index: ', i)) %>%
      colorbar(len = 1)
    
    current_v = default_v
    current_v[5+i] = T
    
    button[[i+1]] = list(method = "restyle",
                         args = list("visible", as.list(current_v)),
                         args2 = list("opacity", set_opacity(data[,i])),
                         label = paste0(i))
  }
  p <- p %>%
    layout(
      updatemenus = list(
        list(buttons = button)))
  
  fig <- subplot(p_lf,p_q1,p_me,p_q3,p_uf, p) 
  fig <- fig %>% layout(title = "Functional boxplot",
                        scene = list(domain=list(x=c(0,0.2),y=c(0.5,1)),
                                     xaxis=axx, yaxis=axx, zaxis=axx,
                                     aspectmode='cube', camera = set_camera),
                        scene2 = list(domain=list(x=c(0.2,0.4),y=c(0.5,1)),
                                      xaxis=axx, yaxis=axx, zaxis=axx,
                                      aspectmode='cube', camera = set_camera),
                        scene3 = list(domain=list(x=c(0.4,0.6),y=c(0.5,1)),
                                      xaxis=axx, yaxis=axx, zaxis=axx,
                                      aspectmode='cube', camera = set_camera),
                        scene4 = list(domain=list(x=c(0.6,0.8),y=c(0.5,1)),
                                      xaxis=axx, yaxis=axx, zaxis=axx,
                                      aspectmode='cube', camera = set_camera), 
                        scene5 = list(domain=list(x=c(0.8,1),y=c(0.5,1)),
                                      xaxis=axx, yaxis=axx, zaxis=axx,
                                      aspectmode='cube', camera = set_camera),
                        scene6 = list(domain=list(x=c(0,1),y=c(0,0.5)),
                                      xaxis=axx, yaxis=axx, zaxis=axx,
                                      aspectmode='cube'),
                        coloraxis=list(colorscale=colorscale, 
                                       cmin = min(c(lf,q3,me,q1,uf)),
                                       cmax = max(c(lf,q3,me,q1,uf)), 
                                       colorbar = list(nticks = 10, tickmode = 'auto',
                                                       len = 0.4, y = .75))
  )
  
  # Print
  print('Click two times on each surface index')
  
  # Display figure
  fig
  
}
outliergram <- function(MEI_data, MBD_data, Fvalue = 1.5, ...){
  # Useful data
  N <- length(MEI_data)
  a_0_2 = -2 / ( N * ( N - 1 ) )
  a_1 = 2 * ( N + 1 ) / ( N - 1 )
  grid_1D = seq( 0, 1, length.out = 100 )
  
  # Compute outliers indices
  d = a_0_2 + a_1 * MEI_data + N^2 * a_0_2 * MEI_data^2 - MBD_data
  Q = stats::quantile( d )
  Q_d3 = Q[ 4 ]
  Q_d1 = Q[ 2 ]
  IQR_d = Q[ 4 ] - Q[ 2 ]
  ID_shape_outlier = which( d >= Q_d3 + Fvalue * IQR_d )
  ID_non_outlying = setdiff(1 : N, ID_shape_outlier)
  
  # Plot 
  fig <- plot_ly() %>%
    add_trace(x = grid_1D, y = a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2, 
              type = 'scatter', mode = 'lines',
              line = list(color = 'darkblue', width = 2, dash = 'dash'), 
              showlegend = F, name = 'Upper parabolic limit') %>%
    add_trace(x = grid_1D, y = a_0_2 + a_1 * grid_1D + a_0_2 * N^2 * grid_1D^2 -
                Q_d3 - Fvalue * IQR_d, type = 'scatter', mode = 'lines',
              line = list(color = 'lightblue', width = 2, dash = 'dash'), 
              showlegend = F, name = 'Lower parabolic limit') %>%
    add_trace(x = MEI_data[ - ID_shape_outlier ], 
              y = MBD_data[ - ID_shape_outlier ], 
              type = 'scatter', mode = 'markers', 
              marker = list(color = "black"), 
              name = "Non outliers") %>%
    add_trace(x = MEI_data[ ID_shape_outlier ], 
              y = MBD_data[ ID_shape_outlier ], 
              type = 'scatter', mode = 'markers', 
              marker = list(color = "red"), 
              name = "Outliers", 
              legendgroup = '1') %>%
    add_text(x = MEI_data[ ID_shape_outlier ],
             y = MBD_data[ ID_shape_outlier ] + 0.5 / 50,
             text = ID_shape_outlier, 
             name = "", 
             legendgroup = '1', showlegend = F) %>%
    layout(xaxis = list(title = 'MEI', range = list(0, 1)), 
           yaxis = list(title = 'MBD', range = list(0, a_0_2 + a_1 / 2 + a_0_2 * N^2/4 )))
  
  fig
  
  return(list(fig,ID_shape_outlier))
}

# Function 
functional_boxplot3D <- function(FEMbasis, lowerfence, q1, median, q3, upperfence, colorscale = 'Set3',
                                   data = NULL, outliergram = FALSE, 
                                   univariate_boxplot = FALSE, persistent = F, ...){
  # Mesh must be 2.5D 
  if(!is(FEMbasis$mesh, "mesh.3D"))
    stop("This function is implemented only for 3D mesh FEM objects")
  
  # Only one plot between outliergram and nodes boxplot profile
  if(outliergram & univariate_boxplot)
    stop("'outliergram' and 'univariate_boxplot' parameters cannot be set both TRUE")
  
  # Must provide data if outliergram is TRUE
  if(outliergram & is.null(data))
    stop("'data' must be provided if outliergram is TRUE")
  
  # Plot 
  p <- plot_functional_boxplot(FEMbasis, lowerfence, q1, median, q3, upperfence, colorscale = colorscale)
  
  ## Functional boxplot
  if(!outliergram & !univariate_boxplot & is.null(data)){
    fig <- p 
  }
  
  ## Functional boxplot with data
  if(!outliergram & !univariate_boxplot & !is.null(data)){
    fig <- plot_functional_boxplot_data(FEMbasis, lowerfence, q1, median, q3, upperfence, data, colorscale = colorscale)
  }
  
  ## Functional boxplot & Outliergram
  if(outliergram){
    # Compute MEI and MBD
    MBD_data <- MBD(t(data))
    MEI_data <- MEI(t(data))
    
    # Set layout options
    set_camera <- list(eye = list(x = 2,y = 2, z = 2))
    axx <- list(
      title = '',
      showgrid = F,
      zeroline = F,
      showticklabels = F 
    )
    
    # Plot functional boxplot and outliergram
    p1 <- outliergram(MEI_data, MBD_data, F_value = 1.5, ...)[[1]]
    p2 <- plot_functional_boxplot_data(FEMbasis, lowerfence, q1, median, q3, upperfence, data, colorscale = colorscale)
    
    # Display 
    fig <- subplot(p2, p1, nrows = 2) %>%
      layout(
        showlegend = F, 
        title = "Functional boxplot & Outliergram",
        scene = list(domain=list(x=c(0,0.2),y=c(0.75,1)),
                     xaxis=axx, yaxis=axx, zaxis=axx,
                     aspectmode='cube', camera = set_camera),
        scene2 = list(domain=list(x=c(0.2,0.4),y=c(0.75,1)),
                      xaxis=axx, yaxis=axx, zaxis=axx,
                      aspectmode='cube', camera = set_camera),
        scene3 = list(domain=list(x=c(0.4,0.6),y=c(0.75,1)),
                      xaxis=axx, yaxis=axx, zaxis=axx,
                      aspectmode='cube', camera = set_camera),
        scene4 = list(domain=list(x=c(0.6,0.8),y=c(0.75,1)),
                      xaxis=axx, yaxis=axx, zaxis=axx,
                      aspectmode='cube', camera = set_camera), 
        scene5 = list(domain=list(x=c(0.8,1),y=c(0.75,1)),
                      xaxis=axx, yaxis=axx, zaxis=axx,
                      aspectmode='cube', camera = set_camera),
        scene6 = list(domain=list(x=c(0,1),y=c(0.5,0.75)),
                      xaxis=axx, yaxis=axx, zaxis=axx,
                      aspectmode='cube')
      )
  }
  
  ## Surface plot & Nodes boxplot profile
  if(univariate_boxplot){
    
    # Create data to be plotted 
    x <- FEMbasis$mesh$nodes[,1] 
    y <- FEMbasis$mesh$nodes[,2] 
    z <- FEMbasis$mesh$nodes[,3] 
    summ <- rbind(lowerfence,q1,median,q3,upperfence)
    N <- ncol(summ)
    df = data.frame()
    for (i in 1:N) {
      temp <- data.frame(rep(x[i],5), rep(y[i],5), rep(z[i],5), summ[,i], rep(i,5))
      df <- rbind(df, temp)
    }
    colnames(df) <- c('x', 'y','z', 'summ', 'id')
    
    # Create shared data
    shared_data <- SharedData$new(df, key = ~id)
    
    # Set selection options 
    s <- attrs_selected(
      opacity = 1,
      marker = list(opacity = 1, color = 'red', size = 8)
    )
    
    # Plot 
    surface <- shared_data %>% 
      plot_ly() %>%
      group_by(id) %>%
      summarise(x = mean(x), y = mean(y), z = mean(z), id = mean(id)) %>%
      add_trace(type = 'mesh3d', 
                x = FEMbasis$mesh$nodes[,1], 
                y = FEMbasis$mesh$nodes[,2], 
                z = FEMbasis$mesh$nodes[,3],
                i = FEMbasis$mesh$triangles[,1]-1,
                j = FEMbasis$mesh$triangles[,2]-1,
                k = FEMbasis$mesh$triangles[,3]-1, 
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'lightgrey',
                opacity = .8,
                hoverinfo = "text",
                text = ~paste("Coordinates: (", round(x,2), round(y,2), round(z,2), ")",
                              "<br>Id: ", id)) %>%
      add_trace(type = 'scatter3d', 
                mode = 'markers',
                x = ~x, y = ~y, z= ~z, size = 1, color = I('black'),
                hoverinfo = "text",
                opacity = 0,
                text = ~paste("")) %>% 
      hide_colorbar() %>%
      layout(scene = list(
        xaxis = list(
          title = '',
          showgrid = F,
          zeroline = F,
          showticklabels = F
        ),
        yaxis = list(
          title = '',
          showgrid = F,
          zeroline = F,
          showticklabels = F
        ),
        zaxis = list(
          title = '',
          showgrid = F,
          zeroline = F,
          showticklabels = F
        ))) 
    
    box <- shared_data %>%
      plot_ly(type = 'box', y = ~summ, color = I('red'), quartilemethod="inclusive", 
              hoverinfo = "y", opacity = 1) %>%
      layout(
        xaxis = list(
          title = '',
          showgrid = F,
          zeroline = F,
          showticklabels = F
        ),
        yaxis = list(
          title = '',
          showgrid = F,
          zeroline = F,
          showticklabels = F
        )
      )
    
    fig <- subplot(surface, box) %>% 
      hide_legend() %>%
      layout(scene = list(domain=list(x=c(0,.5),y=c(0,1)))) %>%
      highlight(on = "plotly_click", off = "plotly_doubleclick", 
                selectize = F, persistent = persistent,
                opacityDim = getOption("opacityDim", 0.5), selected = s)
  }
  
  ## Display figure
  fig
  
}

# Data 

# 1. Hub2.5D built-in data
data("sphere3Ddata")
nodes = sphere3Ddata$nodes
tetrahedrons = sphere3Ddata$tetrahedrons
mesh = create.mesh.3D(nodes = nodes, tetrahedrons = tetrahedrons)
FEMbasis <- create.FEM.basis(mesh)
coeff <- rowSums(nodes)
FEMfunction  <- FEM(coeff, FEMbasis)

# 2. Simulate summary statistics 
beta <- abs(rnorm(1))

me <- beta*coeff
q1 <- beta*coeff - 1 
q3 <- beta*coeff + 1 
lf <- beta*coeff - 2 
uf <- beta*coeff + 3 

# 3. Simulate data on hub nodes

N_good <- 95
N_extra <- 5
N <- N_good + N_extra
data_good <- replicate(n = N_good, expr = rnorm(1)*coeff + rnorm(1,0,3)) 
data_extra <- replicate(n = N_extra, expr = rnorm(1,5,2)*coeff + rnorm(1,0,3)) 
data <- cbind(data_good, data_extra)

# Tests 

# 1. Functional Boxplot
functional_boxplot3D(FEMbasis, lf, q1, me, q3, uf, colorscale = 'Earth')

# 2. Functional Boxplot with data
functional_boxplot3D(FEMbasis, lf, q1, me, q3, uf, data = data, colorscale = 'Earth')

# 3. Functional Boxplot & Outliergram
functional_boxplot3D(FEMbasis, lf, q1, me, q3, uf, data = data, outliergram = T)
