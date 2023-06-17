## FUNCTIONAL BOXPLOT 2D ## 

# Libraries
library(fdaPDE)
library(plotly)
library(crosstalk)
library(roahd)

options(warn = -1)

# Useful functions
plot_functional_boxplot <- function(FEMbasis, f_median, f_q1, f_q3, env_max, env_min){
  
  boundary_ind <- FEMbasis$mesh$segments[,1]
  
  # Create the 'box mesh'
  i <- c(1,sort(rep(2:length(boundary_ind),2)),1)
  j <- i + 1
  j[seq(2, length(j), by = 2)] <- j[seq(2, length(j), by = 2)] + length(j)/2 - 1
  j[length(j)-1] <- 1
  k <- sort(rep((length(j)/2 + 1):length(j),2))
  
  box_triangles = cbind(i,j,k)
  
  plot_ly(type = 'mesh3d',
          x = FEMbasis$mesh$nodes[,1],
          y = FEMbasis$mesh$nodes[,2],
          z = f_median,
          i = FEMbasis$mesh$triangles[,1]-1,
          j = FEMbasis$mesh$triangles[,2]-1,
          k = FEMbasis$mesh$triangles[,3]-1,
          facecolor = rep('grey', length(FEMbasis$mesh$triangles[,1])), 
          showlegend = TRUE, 
          name = 'Median') %>%
    add_trace(type = 'mesh3d',
              x = FEMbasis$mesh$nodes[,1],
              y = FEMbasis$mesh$nodes[,2],
              z = f_q1,
              i = FEMbasis$mesh$triangles[,1]-1,
              j = FEMbasis$mesh$triangles[,2]-1,
              k = FEMbasis$mesh$triangles[,3]-1,
              facecolor = rep('darkgrey', length(FEMbasis$mesh$triangles[,1])), 
              showlegend = TRUE, 
              name = 'Q1') %>%
    add_trace(type = 'mesh3d',
              x = FEMbasis$mesh$nodes[,1],
              y = FEMbasis$mesh$nodes[,2],
              z = f_q3,
              i = FEMbasis$mesh$triangles[,1]-1,
              j = FEMbasis$mesh$triangles[,2]-1,
              k = FEMbasis$mesh$triangles[,3]-1,
              facecolor = rep('darkgrey', length(FEMbasis$mesh$triangles[,1])), 
              showlegend = TRUE, 
              name = 'Q3') %>%
    add_trace(type = 'mesh3d',
              x = FEMbasis$mesh$nodes[,1],
              y = FEMbasis$mesh$nodes[,2],
              z = env_max,
              i = FEMbasis$mesh$triangles[,1]-1,
              j = FEMbasis$mesh$triangles[,2]-1,
              k = FEMbasis$mesh$triangles[,3]-1,
              facecolor = rep('lightgrey', length(FEMbasis$mesh$triangles[,1])), 
              showlegend = TRUE, 
              name = 'Maximum envelope') %>%
    add_trace(type = 'mesh3d',
              x = FEMbasis$mesh$nodes[,1],
              y = FEMbasis$mesh$nodes[,2],
              z = env_min,
              i = FEMbasis$mesh$triangles[,1]-1,
              j = FEMbasis$mesh$triangles[,2]-1,
              k = FEMbasis$mesh$triangles[,3]-1,
              facecolor = rep('lightgrey', length(FEMbasis$mesh$triangles[,1])), 
              showlegend = TRUE, 
              name = 'Minimum envelope') %>%
    add_trace(type = 'mesh3d', 
              x = rep(FEMbasis$mesh$nodes[boundary_ind,1], 2), 
              y = rep(FEMbasis$mesh$nodes[boundary_ind,2], 2), 
              z = c(f_q1[boundary_ind], f_q3[boundary_ind]), 
              i = box_triangles[,1]-1, 
              j = box_triangles[,2]-1, 
              k = box_triangles[,3]-1, 
              opacity = 0.25, 
              name = 'Box') %>%
    layout(
      legend = list(y = 0.5)
    )
  
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
    add_trace(x = MEI_data[ ID_non_outlying ], 
              y = MBD_data[ ID_non_outlying ], 
              type = 'scatter', mode = 'markers', 
              marker = list(color = "black"), 
              hoverinfo = 'text', 
              text = paste("Coordinates:", round(MEI_data[ - ID_shape_outlier ],2), 
                           round(MBD_data[ - ID_shape_outlier ],2), 
                           '<br>Id:', ID_non_outlying),
              name = "Non outliers") %>%
    add_trace(x = MEI_data[ ID_shape_outlier ], 
              y = MBD_data[ ID_shape_outlier ], 
              type = 'scatter', mode = 'markers', 
              marker = list(color = "red"), 
              hoverinfo = 'text', 
              text = paste("Coordinates:", round(MEI_data[ ID_shape_outlier ],2), 
                           round(MBD_data[ ID_shape_outlier ],2), 
                           '<br>Id:', ID_shape_outlier),
              name = "Outliers", 
              legendgroup = '1') %>%
    add_text(x = MEI_data[ ID_shape_outlier ],
             y = MBD_data[ ID_shape_outlier ] + 0.5 / 100,
             text = ID_shape_outlier, 
             name = "", 
             legendgroup = '1', showlegend = F) %>%
    layout(xaxis = list(title = 'MEI', range = list(0, 1)), 
           yaxis = list(title = 'MBD', range = list(0, a_0_2 + a_1 / 2 + a_0_2 * N^2/4 )))
  
  fig
  
  return(list(fig,ID_shape_outlier))
}

# Function
functional_boxplot2D_new <- function(FEMbasis, lowerfence, q1, median, q3, upperfence, 
                                     data = NULL, outliergram = FALSE, univariate_boxplot = FALSE, 
                                     persistent = F, ...){
  # Mesh must be 2D 
  if(!is(FEMbasis$mesh, "mesh.2D"))
    stop("This function is implemented only for 2D mesh FEM objects")
  
  # Only one plot between outliergram and nodes boxplot profile
  if(outliergram & univariate_boxplot)
    stop("'outliergram' and 'univariate_boxplot' parameters cannot be set both TRUE")
  
  # Must provide data if outliergram is TRUE
  if(outliergram & is.null(data))
    stop("'data' must be provided if outliergram is TRUE")
  
  # Plot 
  p <- plot_functional_boxplot(FEMbasis, median, q1, q3, upperfence, lowerfence)
  
  ## Functional boxplot
  if(!outliergram & !univariate_boxplot){
    fig <- p %>%
      layout(title = 'Functional Boxplot')
  }
  
  ## Functional boxplot with outliergram
  if(outliergram){
    
    ### Calcolo outliergram data
    MEI_data <- MEI(data)
    MBD_data <- MBD(data)
    ID_shape_outliers <- outliergram(MEI_data, MBD_data)[[2]]
    
    n_out <- length(ID_shape_outliers)
    
    # Creo il bottone di default in cui è visibile solo il boxplot funzionale
    default <- c(rep(T,6),rep(F,n_out))
    button = list()
    button[[1]] = list(method = "restyle",
                       args = list("visible", as.list(default)),
                       label = "None")
    
    # Update la lista di bottoni e aggiunge la trace corrispondente alla colonna i-esima di funs 
    for (i in 1:length(ID_shape_outliers)) {
      p <- p %>%
        add_trace(type = 'mesh3d',
                  x = FEMbasis$mesh$nodes[,1],
                  y = FEMbasis$mesh$nodes[,2],
                  z = data[ID_shape_outliers[i],],
                  i = FEMbasis$mesh$triangles[,1]-1,
                  j = FEMbasis$mesh$triangles[,2]-1,
                  k = FEMbasis$mesh$triangles[,3]-1,
                  facecolor = rep('red', length(FEMbasis$mesh$triangles[,1])),
                  name = ID_shape_outliers[i], 
                  showlegend = F, 
                  visible = F
        )
      
      current = default
      current[6+i] = T
      
      button[[i+1]] = list(method = "restyle",
                           args = list("visible", as.list(current)),
                           label = paste0(ID_shape_outliers[i]))
      
      
    }
    
    p <- p %>%
      layout(
        updatemenus = list(
          list(buttons = button)))
    
    fig <- subplot(p,outliergram(MEI_data, MBD_data)[[1]])
    fig <- fig %>%
      layout(
        scene = list(
          domain = list(x = c(0,0.5), y = c(0,1))
        ), 
        scene2 = list(
          domain = list(x = c(0.5,1), y = c(0,1))
        ), 
        showlegend = FALSE,
        showlegend2 = FALSE)
  }
  
  ## Profilo boxplot univariato per ogni nodo 
  if(univariate_boxplot){
    # Layout options 
    axx <- list(
      title = '',
      showgrid = F,
      zeroline = F,
      showticklabels = F
    )
    
    # Create data to be plotted 
    x <- FEMbasis$mesh$nodes[,1] 
    y <- FEMbasis$mesh$nodes[,2]  
    summ <- rbind(lowerfence,q1,median,q3,upperfence)
    N <- ncol(summ)
    df = data.frame()
    for (i in 1:N) {
      temp <- data.frame(rep(x[i],5), rep(y[i],5), summ[,i], rep(i,5))
      df <- rbind(df, temp)
    }
    colnames(df) <- c('x', 'y', 'summ', 'id')
    
    # Create shared data
    shared_data <- SharedData$new(df, key = ~id)
    
    # Set selection options 
    s <- attrs_selected(
      marker = list(opacity = 1, color = 'red', size = 8)
    )
    
    # Plot 
    
    
    
    box <- shared_data %>%
      plot_ly(type = 'box', y = ~summ, color = I('red'), quartilemethod="inclusive", 
              hoverinfo = "y", opacity = 1) %>%
      layout(xaxis = list(title = '',
                          zerolinecolor = '#ffff',
                          zerolinewidth = 2,
                          gridcolor = 'ffff',
                          showticklabels=FALSE), 
             yaxis = list(title = '',
                          zerolinecolor = '#ffff',
                          zerolinewidth = 2,
                          gridcolor = 'ffff',
                          showticklabels=FALSE)) 
    
    boundary_ind <- FEMbasis$mesh$segments[,1]
    
    # Create the 'box mesh'
    i <- c(1,sort(rep(2:length(boundary_ind),2)),1)
    j <- i + 1
    j[seq(2, length(j), by = 2)] <- j[seq(2, length(j), by = 2)] + length(j)/2 - 1
    j[length(j)-1] <- 1
    k <- sort(rep((length(j)/2 + 1):length(j),2))
    
    box_triangles = cbind(i,j,k)
    
    func_boxplot <- shared_data %>% 
      plot_ly(scene = 'scene2') %>%
      group_by(id) %>%
      summarise(x = mean(x), y = mean(y), id = mean(id)) %>%
      add_trace(type = 'mesh3d', 
                x = FEMbasis$mesh$nodes[,1], 
                y = FEMbasis$mesh$nodes[,2], 
                z = me,
                i = FEMbasis$mesh$triangles[,1]-1,
                j = FEMbasis$mesh$triangles[,2]-1,
                k = FEMbasis$mesh$triangles[,3]-1, 
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'grey',
                opacity = 1,
                hoverinfo = "text",
                text = ~paste("Coordinates: (", round(x,2), round(y,2), round(me,2), ")",
                              "<br>Id: ", id)) %>%
      add_trace(type = 'scatter3d', 
                mode = 'markers',
                x = ~x, y = ~y,z = me, size = .01, opacity = 0,
                hoverinfo = "text",
                text = ~paste("")) %>% 
      add_trace(type = 'mesh3d', 
                x = FEMbasis$mesh$nodes[,1], 
                y = FEMbasis$mesh$nodes[,2], 
                z = lf,
                i = FEMbasis$mesh$triangles[,1]-1,
                j = FEMbasis$mesh$triangles[,2]-1,
                k = FEMbasis$mesh$triangles[,3]-1, 
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'grey',
                opacity = 1,
                hoverinfo = "text",
                text = ~paste("Coordinates: (", round(x,2), round(y,2), round(lf,2), ")",
                              "<br>Id: ", id)) %>%
      add_trace(type = 'scatter3d', 
                mode = 'markers',
                x = ~x, y = ~y,z = lf, size = .01, opacity = 0,
                hoverinfo = "text",
                text = ~paste("")) %>% 
      add_trace(type = 'mesh3d', 
                x = FEMbasis$mesh$nodes[,1], 
                y = FEMbasis$mesh$nodes[,2], 
                z = uf,
                i = FEMbasis$mesh$triangles[,1]-1,
                j = FEMbasis$mesh$triangles[,2]-1,
                k = FEMbasis$mesh$triangles[,3]-1, 
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'lightgrey',
                opacity = 1,
                hoverinfo = "text",
                text = ~paste("Coordinates: (", round(x,2), round(y,2), round(uf,2), ")",
                              "<br>Id: ", id)) %>%
      add_trace(type = 'scatter3d', 
                mode = 'markers',
                x = ~x, y = ~y,z = uf, size = .01, opacity = 0,
                hoverinfo = "text",
                text = ~paste("")) %>% 
      add_trace(type = 'mesh3d', 
                x = FEMbasis$mesh$nodes[,1], 
                y = FEMbasis$mesh$nodes[,2], 
                z = q1,
                i = FEMbasis$mesh$triangles[,1]-1,
                j = FEMbasis$mesh$triangles[,2]-1,
                k = FEMbasis$mesh$triangles[,3]-1, 
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'darkgrey',
                opacity = .6,
                hoverinfo = "text",
                text = ~paste("Coordinates: (", round(x,2), round(y,2), round(q1,2), ")",
                              "<br>Id: ", id)) %>%
      add_trace(type = 'scatter3d', 
                mode = 'markers',
                x = ~x, y = ~y,z = q1, size = .01, opacity = 0,
                hoverinfo = "text",
                text = ~paste("")) %>% 
      add_trace(type = 'mesh3d', 
                x = FEMbasis$mesh$nodes[,1], 
                y = FEMbasis$mesh$nodes[,2], 
                z = q3,
                i = FEMbasis$mesh$triangles[,1]-1,
                j = FEMbasis$mesh$triangles[,2]-1,
                k = FEMbasis$mesh$triangles[,3]-1, 
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'darkgrey',
                opacity = .6,
                hoverinfo = "text",
                text = ~paste("Coordinates: (", round(x,2), round(y,2), round(q3,2), ")",
                              "<br>Id: ", id)) %>%
      add_trace(type = 'scatter3d', 
                mode = 'markers',
                x = ~x, y = ~y,z = q3, size = .01, opacity = 0,
                hoverinfo = "text",
                text = ~paste("")) %>% 
      add_trace(type = 'mesh3d', 
                x = rep(FEMbasis$mesh$nodes[boundary_ind,1], 2), 
                y = rep(FEMbasis$mesh$nodes[boundary_ind,2], 2), 
                z = c(q1[boundary_ind], q3[boundary_ind]), 
                i = box_triangles[,1]-1, 
                j = box_triangles[,2]-1, 
                k = box_triangles[,3]-1,
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'grey',
                opacity = 0.25, 
                name = 'Box') %>%
      add_trace(type = 'mesh3d', 
                x = FEMbasis$mesh$nodes[,1], 
                y = FEMbasis$mesh$nodes[,2], 
                z = rep(min(lf) - (max(uf) - min(lf))/2,N),
                i = FEMbasis$mesh$triangles[,1]-1,
                j = FEMbasis$mesh$triangles[,2]-1,
                k = FEMbasis$mesh$triangles[,3]-1, 
                intensity = rep(0,N),
                color = rep(0,N),
                colors = 'lightgrey',
                opacity = .8,
                hoverinfo = "text",
                text = ~paste("Coordinates: (", round(x,2), round(y,2), ")",
                              "<br>Id: ", id)) %>%
      add_trace(type = 'scatter3d', 
                mode = 'markers',
                x = ~x, y = ~y,z = rep(min(lf) - (max(uf) - min(lf))/2,N), 
                size = .1,
                opacity = 0, 
                color = I('grey'),
                hoverinfo = "text",
                text = ~paste("")) %>% 
      hide_colorbar() %>%
      hide_legend() %>%
      layout(scene = list(
        aspectmode = 'data',
        xaxis = axx,
        yaxis = axx,
        zaxis = axx)) 
    
    
    
    fig <- subplot(func_boxplot, box) %>% 
      hide_legend() %>%
      layout(scene2 = list(domain=list(x=c(0,.5),y=c(0,1)), 
                           xaxis = axx,
                           yaxis = axx,
                           zaxis = axx), 
             xaxis = list(domain=list(x=c(0,0.5),y=c(0,.5)))) %>%
      highlight(on = "plotly_click", off = "plotly_doubleclick", 
                selectize = F, persistent = F,
                opacityDim = getOption("opacityDim", 0.8), selected = s)
    
    
    
  }
  
  # Display the figure
  fig
  
}

# Simulated Data 

# 1. Horseshoe2D built-in data
data(horseshoe2D)
boundary_nodes = horseshoe2D$boundary_nodes
boundary_segments = horseshoe2D$boundary_segments
locations = horseshoe2D$locations
mesh = create.mesh.2D(nodes = rbind(boundary_nodes, locations), segments = boundary_segments)
FEMbasis = create.FEM.basis(mesh)
coeff = fs.test(mesh$nodes[,1], mesh$nodes[,2]) 
FEMfunction = FEM(coeff, FEMbasis)

# 2. Simulate data on horseshoe2D domain
N_good <- 290
N_extra <- 10
N <- N_good + N_extra
data_good <- replicate(n = N_good, expr = rnorm(1)*coeff + rnorm(1,0,3)) 
id_good <- 1:290
data_extra <- replicate(n = N_extra, expr = rnorm(1,5,1)*coeff + rnorm(1,0,3))
id_extra <- 291:300
data <- t(cbind(data_good, data_extra))

# 3. Simulate summary statistics
me = coeff 
q1 = coeff - 1
q3 = coeff + 1
uf = coeff + 2
lf = coeff - 2

# Tests 

# Test 1: Functional Boxplot
functional_boxplot2D_new(FEMbasis, lf, q1, me, q3, uf) 

# Test 2: Functional Boxplot & Outliergram
functional_boxplot2D_new(FEMbasis, lf, q1, me, q3, uf, data = data, outliergram = T) 

# Test 3: Image plot & Univariate boxplot
functional_boxplot2D_new(FEMbasis, lf, q1, me, q3, uf, univariate_boxplot = T)



