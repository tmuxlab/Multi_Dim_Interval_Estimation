#点推定と信頼区間
loc_grid <- seq(0, 1, length.out = 256)

mat_hat <- vec1024_to_mat256x4(bs$lambdaZ_hat_T)
mat_lo  <- vec1024_to_mat256x4(bs$lambdaZ_CI_T$lo)
mat_hi  <- vec1024_to_mat256x4(bs$lambdaZ_CI_T$hi)

par(mfrow=c(2,2), mar=c(4,4,2,1))
for(mod in 1:4){
        yl <- range(mat_lo[,mod], mat_hi[,mod], finite = TRUE)
        plot(loc_grid, mat_hat[,mod], type="l",ylim = yl,
             xlab="LOCnorm", ylab=expression(lambda[z](T,z)),
             main=paste0("module ", mod-1))
        lines(loc_grid, mat_lo[,mod], lty=2)
        lines(loc_grid, mat_hi[,mod], lty=2)
}
par(mfrow=c(1,1))

#マーク断面（m=6.0）
k <- which(bs$m_set == 6.0)
mat_hat_m <- vec1024_to_mat256x4(bs$lambdaZM_hat_T[[k]])
mat_lo_m  <- vec1024_to_mat256x4(bs$lambdaZM_CI_T[[k]]$lo)
mat_hi_m  <- vec1024_to_mat256x4(bs$lambdaZM_CI_T[[k]]$hi)

par(mfrow=c(2,2), mar=c(4,4,2,1))
for(mod in 1:4){
        yl <- range(mat_lo_m[,mod], mat_hi_m[,mod], finite = TRUE)
        plot(loc_grid, mat_hat_m[,mod], type="l",ylim = yl,
             xlab="LOCnorm", ylab=expression(lambda(T,z,m)),
             main=paste0("module ", mod-1, " (m=", bs$m_set[k], ")"))
        lines(loc_grid, mat_lo_m[,mod], lty=2)
        lines(loc_grid, mat_hi_m[,mod], lty=2)
}
par(mfrow=c(1,1))
#####
par(mfrow=c(2,2), mar=c(4,4,2,1))

for(mod in 1:4){
        
        yvals <- c(mat_hat_m[,mod], mat_lo_m[,mod], mat_hi_m[,mod])
        yvals <- yvals[is.finite(yvals)]
        
        yl <- range(yvals)                 # 自動レンジ
        yt <- pretty(yl, n = 5)            # 見やすい目盛りを自動生成（だいたい5本）
        
        plot(loc_grid, mat_hat_m[,mod],
             type="l", ylim=range(yt), yaxt="n",
             xlab="LOCnorm", ylab=expression(lambda(T,z,m)),
             main=paste0("module ", mod-1, " (m=", bs$m_set[k], ")"))
        
        axis(2, at=yt, labels=format(yt, digits=3))  # 桁数も適度に
        lines(loc_grid, mat_lo_m[,mod], lty=2)
        lines(loc_grid, mat_hi_m[,mod], lty=2)
}

par(mfrow=c(1,1))
######
#透明区間
par(mfrow=c(2,2), mar=c(4,4,2,1))
for(mod in 1:4){
        lo <- mat_lo[,mod]; hi <- mat_hi[,mod]; hat <- mat_hat[,mod]
        okp <- is.finite(loc_grid) & is.finite(lo) & is.finite(hi) & is.finite(hat)
        
        yl <- range(lo[okp], hi[okp])
        plot(loc_grid[okp], hat[okp], type="n", ylim=yl,
             xlab="LOCnorm", ylab=expression(lambda[z](T,z)),
             main=paste0("module ", mod-1))
        
        polygon(c(loc_grid[okp], rev(loc_grid[okp])),
                c(lo[okp], rev(hi[okp])),
                col = rgb(0,0,0,0.12), 
                border = NA)
        
        lines(loc_grid[okp], hat[okp], lwd=2)
}
par(mfrow=c(1,1))
#網掛け
# par(mfrow=c(2,2), mar=c(4,4,2,1))
# for(mod in 1:4){
#         lo <- mat_lo[,mod]; hi <- mat_hi[,mod]; hat <- mat_hat[,mod]
#         okp <- is.finite(loc_grid) & is.finite(lo) & is.finite(hi) & is.finite(hat)
#         
#         yl <- range(lo[okp], hi[okp])
#         plot(loc_grid[okp], hat[okp], type="n", ylim=yl,
#              xlab="LOCnorm", ylab=expression(lambda[z](T,z)),
#              main=paste0("module ", mod-1))
#         
#         polygon(c(loc_grid[okp], rev(loc_grid[okp])),
#                 c(lo[okp], rev(hi[okp])),
#                 density = 20, angle = 45,  # ★網掛け
#                 border = NA)
#         
#         lines(loc_grid[okp], hat[okp], lwd=2)
#         lines(loc_grid[okp], lo[okp], lty=2)
#         lines(loc_grid[okp], hi[okp], lty=2)
# }
# par(mfrow=c(1,1))
