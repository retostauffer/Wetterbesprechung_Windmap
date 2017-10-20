# -------------------------------------------------------------------
# - NAME:        georg.R
# - AUTHOR:      Reto Stauffer
# - DATE:        2017-10-20
# -------------------------------------------------------------------
# - DESCRIPTION:
# -------------------------------------------------------------------
# - EDITORIAL:   2017-10-20, RS: Created file on thinkreto.
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2017-10-20 19:14 on thinkreto
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Helper function to convert uv to dd and ff
# -------------------------------------------------------------------
   uv2ddff <- function( u, v, rad = F ) {
     ## polar coordinates:
     ff <- sqrt(u^2 + v^2)
     dd <- atan(v/u) + (u < 0) * pi
     idx <- which( dd < 0 );    dd[ idx ] <- dd[ idx ] + 2 * pi
     ## convert angle to meteorological convention
     dd <- 3 * pi / 2 - dd
     idx <- which( dd < 0 );    dd[ idx ] <- dd[ idx ] + 2 * pi
     ## if rad (radiants) = F we have to convert to degrees.
     if ( ! rad ) dd <- dd * 180 / pi
     data.frame(ff, dd)
   }


# -------------------------------------------------------------------
# Draw longitude latitude grid
# X: proj4string of 
# -------------------------------------------------------------------
   drawLonLatGrid <- function( newproj, lats = seq(0,85,by=5), lons = c(-30,30,by=5), n=360, ... ) {
      stopifnot( is.character( newproj ) )
      llProj4string  <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"
      # Grid lines?
      for ( lat in seq(0,85,by=5) ) {
         ll <- SpatialPoints( data.frame( x = seq(-180,180,length=n), y = rep(lat,n) ),
               proj4string = CRS(llProj4string) )
         ll <- spTransform( ll, crs(data) )
         lines( ll@coords[,1], ll@coords[,2], ... )
      }
      for ( lon in seq(-30,30,by=5) ) {
         ll <- SpatialPoints( data.frame( x = rep(lon, n), y = seq(-90,90,length=n) ),
               proj4string = CRS(llProj4string) )
         ll <- spTransform( ll, crs(data) )
         lines( ll@coords[,1], ll@coords[,2], ... )
      }
   }

# -------------------------------------------------------------------
# Adding somehow meaningful axis labels
# data: the raster object (we only need the proj4string out of it)
# -------------------------------------------------------------------
   addAxis <- function( data ) {
      llProj4string  <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"
      # Try to add useful axis
      xax <- pretty(as.matrix(data@extent)["x",])
      sp  <- SpatialPoints( data.frame("x"=xax,"y"=rep(attr(data@extent,"ymin"),length(xax))),
             proj4string = crs( data ) )
      sp  <- spTransform( sp, llProj4string )
      axis( side = 1, at = xax, labels = sprintf("%.1f",sp@coords[,"x"] ) )
   
      yax <- pretty(as.matrix(data@extent)["y",])
      sp  <- SpatialPoints( data.frame("x"=yax,"y"=rep(attr(data@extent,"xmin"),length(yax))),
             proj4string = crs( data ) )
      sp  <- spTransform( sp, llProj4string )
      axis( side = 2, at = yax, labels = sprintf("%.1f",sp@coords[,"y"] ) )
   }

# -------------------------------------------------------------------
# Helper function to load grib data as raster and reproject the
# data. Note that I am only reading grib message 1, wherefore
# I've subsetted the grib files such that I know what's on message 1.
# -------------------------------------------------------------------
   getrasterdata <- function( x, proj4string = NULL ) {
      stopifnot( file.exists(x) )
      require("raster")
      r <- raster( x )
      if ( is.null(new) ) return(r)
      r <- projectRaster( r, crs = CRS(proj4string) )
   
   }


# -------------------------------------------------------------------
# Pretty hardcoded function which loads a set of grib files and
# returns one named RasterStack object with u/v/z/dd/ff.
# -------------------------------------------------------------------
   prepareRasterData <- function( newProj4string, extent = NULL ) {
      if ( missing(newProj4string) )
         newProj4string <- "+proj=stere +lat_0=90 +lon_0=-10 +k=1 +x_0=0 +y_0=-0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
      u500 <- getrasterdata( "demodata/ECMWF_012_u500.grib", newProj4string )
      v500 <- getrasterdata( "demodata/ECMWF_012_v500.grib", newProj4string )
      z500 <- getrasterdata( "demodata/ECMWF_012_z500.grib", newProj4string )
      
      dd500 <- ff500 <- u500 # Clone raster
      values(ff500) <- uv2ddff( values(u500), values(v500) )$ff
      values(dd500) <- uv2ddff( values(u500), values(v500) )$dd
      
      data <- stack( u500, v500, z500, ff500, dd500 )
      names(data) <- c("u500","v500","z500","ff500","dd500")
      if ( ! is.null(extent) ) data <- crop( data, extent )
      data
   }


# -------------------------------------------------------------------
# Function to create spatial line barbs 
# -------------------------------------------------------------------
   gentemplate_barb_template <- function( x ) {
      x <- round( x / 5 ) * 5
   
      add50 <- function( x ) {
         y <- 1 - (nrow(x)-1)*0.1
         x <- rbind( x, c( 0.0, y,     0.4, y-0.1) )
         x <- rbind( x, c( 0.0, y-0.2, 0.4, y-0.1) )
      }
      add10 <- function( x ) {
         y <- 1 - (nrow(x)-1)*0.1
         x <- rbind( x, c( 0.0, y, 0.4, y+0.1) )
      }
      add5 <- function( x ) {
         y <- 1 - (nrow(x)-1)*0.1
         x <- rbind( x, c( 0.0, y, 0.2, y+0.05) )
      }
   
      res <- matrix( c(0.0, 0.0, 0.0, 1.0), ncol = 4, dimnames = list(NULL,c("x0","y0","x1","y1")) )
      while ( x > 0 ) {
         if ( floor(x/50) > 0 ) {
            res <- add50(res);    x <- x - 50
         } else if ( floor(x/10) > 0 ) {
            res <- add10(res);    x <- x - 10
         } else if ( x == 5 ) {
            res <- add5(res);     x <- x - 5
         }
      }
   
      matrix( as.numeric( t(res) ), ncol = 2, byrow = T,
              dimnames = list(NULL,c("x","y")))
   }
   
   
   generate_barbs <- function( x, y, dd, ff, scale ) {

      w <- par("pin")[1]/diff(par("usr")[1:2])
      h <- par("pin")[2]/diff(par("usr")[3:4])
      asp <- w/h
   
      # Meteorological wind direction in math radiant
      fun <- function(dd) (360-dd) / 180 * pi

      # Only consider the ones where dd is not equal to NA
      res <- list()
      idx <- which(! is.na(dd) & ! is.na(ff) )
      for ( i in idx ) {
   
         # Wind direction in radiant
         ddrad <- fun(dd[i])
   
         # Rotation matrix
         R <- matrix( c(cos(ddrad),-sin(ddrad),sin(ddrad),cos(ddrad)), ncol = 2)
         R[,1] <- R[,1] * asp # Aspect for elliptical rotation
         R[,2] <- R[,2]
   
         # Scale, rotate and re-center to x/y
         xx     <- (scale * gentemplate_barb_template( ff[i] )) %*% R
         xx[,1] <- xx[,1] + x[i]
         xx[,2] <- xx[,2] + y[i]
   
         # Create polygon
         p <- list()
         for ( k in 1:(nrow(xx)/2) ) {
            p[[k]] <- sp::Line( xx[ (2*k)-c(1,0), ] )
         }
         p <- sp::Lines( p, i )
         # Append
         res[[length(res)+1]]   <- p
   
      }
      sp::SpatialLines( res )
   }
   

# -------------------------------------------------------------------
# MAIN SCRIPT : MAIN SCRIPT : MAIN SCRIPT : MAIN SCRIPT : MAIN SCRIPT
# -------------------------------------------------------------------

   library("sp"); library("raster")

   llProj4string  <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +no_defs +towgs84=0,0,0"
   newProj4string <- "+proj=stere +lat_0=90 +lon_0=10 +k=1 +x_0=0 +y_0=-0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
   # Define two points (long/lat) which define the crop area
   e <- SpatialPoints( data.frame( x = c(-10,20), y = c(30,50) ),
         proj4string = CRS( llProj4string ) )
   e <- spTransform( e, newProj4string )
   e <- extent( t(e@coords) )

   data <- prepareRasterData( newProj4string, e )

   # In addition we need a world map. Use Polygon data set.
   # Crop countries with the same extent (take it from data@extent, is e)
   library('rworldxtra')
   data(countriesHigh)
   countriesHigh <- crop( spTransform( countriesHigh, newProj4string ), data@extent )


# -------------------------------------------------------------------
# This is the automatic version
# -------------------------------------------------------------------

   # Plotting base map
   plot( data$ff500, main = "Georg's Dingsda", xaxt = "n", yaxt = "n" )
   plot( countriesHigh, border = "gray50", add = TRUE )
   contour( data$z500/9.806, add = TRUE )

   # Draw longitude latitude grid
   drawLonLatGrid( proj4string(data), col = 4 )

   # Adding somehow meaningful axis labels
   addAxis( data )

   # Picking random coordinates to draw the barbs
   idx <- sample(1:prod(dim(data)[1:2]),200)
   barbs <- list( x   = coordinates(data)[idx,1],
                  y   = coordinates(data)[idx,2],
                  dd  = values(data$dd500)[idx],
                  ff  = values(data$ff500)[idx] * 1.94384 )

   # Draw barbs
   spbarbs <- generate_barbs( barbs$x, barbs$y, barbs$dd, barbs$ff, 1e5 )
   plot( spbarbs, add = TRUE, lwd = 2 )
   points( barbs$x, barbs$y, pch = 19, cex = 0.5 )



# -------------------------------------------------------------------
# This is the automatic version
# -------------------------------------------------------------------

   stop(" from here on - manual version (manually select barbs " )

   # Plotting base map
   plot( data$ff500, main = "Georg's Dingsda", xaxt = "n", yaxt = "n" )
   plot( countriesHigh, border = "gray50", add = TRUE )
   contour( data$z500/9.806, add = TRUE )
   drawLonLatGrid( proj4string(data), col = 4 )
   addAxis( data )

   # Manual select
   n <- 10
   barbs <- list()
   while ( length(barbs$x) < n ) {
      cat(sprintf(" Please select barb number %d/%d: ", length(barbs$x)+1,n))
      tmp.xy   <- locator(1)
      tmp.xy   <- SpatialPoints(data.frame(x=tmp.xy$x,y=tmp.xy$y),proj4string=crs(data))
      tmp.data <- as.data.frame( extract( data, tmp.xy ) )
      if ( any(is.na(tmp.data)) ) cat("Sorry, out of defined area ...\n")
      # Else pick x/y
      barbs$x  <- c(barbs$x,  tmp.xy$x)
      barbs$y  <- c(barbs$y,  tmp.xy$y)
      barbs$dd <- c(barbs$dd, tmp.data$dd500)
      barbs$ff <- c(barbs$ff, tmp.data$ff500)
      cat("Well done! Next please!\n")
   }

   # Draw barbs
   spbarbs <- generate_barbs( barbs$x, barbs$y, barbs$dd, barbs$ff, 1e5 )
   plot( spbarbs, add = TRUE, lwd = 2 )
   points( barbs$x, barbs$y, pch = 19, cex = 0.5 )















