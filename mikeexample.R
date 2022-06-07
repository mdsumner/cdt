library(silicate)
sc <- SC0(inlandwaters)

#f <- raadfiles::thelist_files(pattern = "parcel.*hobart")
#sc <- silicate::SC0(sf::read_sf(f$fullname))

library(cdt)  
del <- head(tridel(as.matrix(sc_vertex(sc)), 
            as.matrix(sc_edge(sc)[c(".vx0", ".vx1")]), mark_domains = F), -1)



#remotes::install_github("hypertidy/plover")

## triangle centroids
tric <- cbind(colMeans(matrix(sc_vertex(sc)$x_[t(del)], 3L)), 
      colMeans(matrix(sc_vertex(sc)$y_[t(del)], 3L)))


idx <- plover::ic_over(tric, inlandwaters$geom)
## bad triangles we chuck
del <- del[!is.na(idx), ]
range(idx)
range(del)
plot(sc_vertex(sc), pch = ".")
#plot(tric, pch = 19, cex = .2)
lines(sc_vertex(sc)[c(t(cbind(del[,c(1:3, 1)], NA))), ], col = "firebrick")
