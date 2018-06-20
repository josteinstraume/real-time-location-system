txt = readLines("../data/offline_final_trace.txt");
# sum(substr(txt, 1, 1) == "#");
# length(txt);
# strsplit(txt[4], ";")[[1]];
# unlist(lapply(strsplit(txt[4], ";")[[1]], function(x) sapply(strsplit(x, "=")[[1]], strsplit, ",")));
# tokens = strsplit(txt[4], "[;=,]")[[1]];
# #tokens[1:10];
# #tokens[c(2, 4, 6:8, 10)]
# 
# 
# 
# tokens[ - ( 1:10 ) ]
# tmp = matrix(tokens[ - (1:10) ], ncol = 4, byrow = TRUE)
# mat = cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp), ncol = 6, byrow = TRUE), tmp)
# dim(mat)
# processLine =
#   function(x)
#   {
#     tokens = strsplit(x, "[;=,]")[[1]]
#     tmp = matrix(tokens[ - (1:10) ], ncol = 4, byrow = TRUE)
#     cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp),
#                  ncol = 6, byrow = TRUE), tmp)
#   }
# tmp = lapply(txt[4:20], processLine)
# sapply(tmp, nrow)
# offline = as.data.frame(do.call("rbind", tmp))
# dim(offline)
# 
lines = txt[ substr(txt, 1, 1) != "#" ]
# options(error = recover, warn = 2)
# tmp = lapply(lines, processLine)


processLine = function(x) {
  tokens = strsplit(x, "[;=,]")[[1]]
  if (length(tokens) == 10)
    return(NULL)
  tmp = matrix(tokens[ - (1:10) ], ncol = 4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow(tmp), 6, byrow = TRUE), tmp)
}

options(error = recover, warn = 1)
tmp = lapply(lines, processLine)
offline = as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)
dim(offline)

names(offline) = c("time", "scanMac", "posX", "posY", "posZ",
                   "orientation", "mac", "signal",
                   "channel", "type")

numVars = c("time", "posX", "posY", "posZ",
            "orientation", "signal")

offline[ numVars ] =  lapply(offline[ numVars ], as.numeric)
offline = offline[ offline$type == "3", ]
offline = offline[ , "type" != names(offline) ]
dim(offline)
offline$rawTime = offline$time
offline$time = offline$time/1000
class(offline$time) = c("POSIXt", "POSIXct")
unlist(lapply(offline, class))
summary(offline[, numVars])
summary(sapply(offline[ , c("mac", "channel", "scanMac")],
               as.factor))
offline = offline[ , !(names(offline) %in% c("scanMac", "posZ"))]
length(unique(offline$orientation))
plot(ecdf(offline$orientation))
roundOrientation = function(angles) {
  refs = seq(0, by = 45, length  = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}
offline$angle = roundOrientation(offline$orientation)
with(offline, boxplot(orientation ~ angle,
                      xlab = "nearest 45 degree angle",
                      ylab="orientation"))
c(length(unique(offline$mac)), length(unique(offline$channel)))
table(offline$mac)
subMacs = names(sort(table(offline$mac), decreasing = TRUE))[1:7]
offline = offline[ offline$mac %in% subMacs, ]
macChannel = with(offline, table(mac, channel))
apply(macChannel, 1, function(x) sum(x > 0))
offline = offline[ , "channel" != names(offline)]
locDF = with(offline,
             by(offline, list(posX, posY), function(x) x))
length(locDF)
sum(sapply(locDF, is.null))
locDF = locDF[ !sapply(locDF, is.null) ]
length(locDF)
locCounts = sapply(locDF, nrow)
locCounts = sapply(locDF,
                   function(df)
                     c(df[1, c("posX", "posY")], count = nrow(df)))
class(locCounts)
dim(locCounts)
locCounts[ , 1:8]

locCounts = t(locCounts)
plot(locCounts, type = "n", xlab = "", ylab = "")
text(locCounts, labels = locCounts[,3], cex = .8, srt = 45)


#offlineRedo = readData()
#identical(offline, offlineRedo)


library(lattice)
bwplot(signal ~ factor(angle) | mac, data = offline,
       subset = posX == 2 & posY == 12
       & mac != "00:0f:a3:39:dd:cd",
       layout = c(2,3))

summary(offline$signal)

densityplot( ~ signal | mac + factor(angle), data = offline,
             subset = posX == 24 & posY == 4 &
               mac != "00:0f:a3:39:dd:cd",
             bw = 0.5, plot.points = FALSE)

offline$posXY = paste(offline$posX, offline$posY, sep = "-")

byLocAngleAP = with(offline,
                    by(offline, list(posXY, angle, mac),
                       function(x) x))

signalSummary = lapply(byLocAngleAP, function(oneLoc)
        {
           ans = oneLoc[1, ]
           ans$medSignal = median(oneLoc$signal)
           ans$avgSignal = mean(oneLoc$signal)
           ans$num = length(oneLoc$signal)
           ans$sdSignal = sd(oneLoc$signal)
           ans$iqrSignal = IQR(oneLoc$signal)
           ans
         })

offlineSummary = do.call("rbind", signalSummary)

breaks = seq(-90, -30, by = 5)
bwplot(sdSignal ~ cut(avgSignal, breaks = breaks),
       data = offlineSummary,
       subset = mac != "00:0f:a3:39:dd:cd",
       xlab = "Mean Signal", ylab = "SD Signal")

with(offlineSummary,
     smoothScatter((avgSignal - medSignal) ~ num,
                   xlab = "Number of Observations",
                   ylab = "mean - median"))
abline(h = 0, col = "#984ea3", lwd = 2)


lo.obj =
  with(offlineSummary,
       loess(diff ~ num,
             data = data.frame(diff = (avgSignal - medSignal),
                               num = num)))

lo.obj.pr = predict(lo.obj, newdata = data.frame(num = (70:120)))
lines(x = 70:120, y = lo.obj.pr, col = "#4daf4a", lwd = 2)

oneAPAngle = subset(offline, mac == subMacs[5] & angle == 0)
oneAPAngle = subset(offlineSummary,
                    mac == subMacs[5] & angle == 0)

install.packages("fields")
library(fields)
smoothSS = Tps(oneAPAngle[, c("posX","posY")],
               oneAPAngle$avgSignal)

vizSmooth = predictSurface(smoothSS)
plot.surface(vizSmooth, type = "C")
points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)

surfaceSS = function(data, mac, angle) {
  #newDat = subset(dat, dat$mac == macX & dat$angle == angleX)
  oneAPAngle = data[ data$mac == mac & data$angle == angle, ]
  smoothSS = Tps(oneAPAngle[, c("posX","posY")],
                 oneAPAngle$avgSignal)
  
  vizSmooth = predictSurface(smoothSS)
  plot.surface(vizSmooth, type = "C")
  points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)
  #return(newDat)
}
parCur = par(mfrow = c(2,2), mar = rep(1, 4))

mapply(surfaceSS, mac = subMacs[ rep(c(5, 1), each = 2) ],
       angle = rep(c(0, 135), 2),
       data = list(data = offlineSummary))

offlineSummary = subset(offlineSummary, mac != subMacs[2])

AP = matrix( c( 7.5, 6.3, 2.5, -.8, 12.8, -2.8,
                1, 14, 33.5, 9.3,  33.5, 2.8),
             ncol = 2, byrow = TRUE,
             dimnames = list(subMacs[ -2 ], c("x", "y") ))

diffs = offlineSummary[ , c("posX", "posY")] -
  AP[ offlineSummary$mac, ]


offlineSummary$dist = sqrt(diffs[ , 1]^2 + diffs[ , 2]^2)

xyplot(signal ~ dist | factor(mac) + factor(angle),
       data = offlineSummary, pch = 19, cex = 0.3,
       xlab ="distance")

oldPar = par(mar = c(3.1, 3.1, 1, 1))
library(lattice)
xyplot(signal ~ dist | factor(mac) + factor(angle), data = offlineSummary, pch = 19, cex = 0.3, xlab ="distance")
par(oldPar)
readData = 
  function(filename = 'Data/offline.final.trace.txt', 
           subMacs = c("00:0f:a3:39:e1:c0", "00:0f:a3:39:dd:cd", "00:14:bf:b1:97:8a",
                       "00:14:bf:3b:c7:c6", "00:14:bf:b1:97:90", "00:14:bf:b1:97:8d",
                       "00:14:bf:b1:97:81"))
  {
    txt = readLines(filename)
    lines = txt[ substr(txt, 1, 1) != "#" ]
    tmp = lapply(lines, processLine)
    offline = as.data.frame(do.call("rbind", tmp), 
                            stringsAsFactors= FALSE) 
    
    names(offline) = c("time", "scanMac", 
                       "posX", "posY", "posZ", "orientation", 
                       "mac", "signal", "channel", "type")
    
    # keep only signals from access points
    offline = offline[ offline$type == "3", ]
    
    # drop scanMac, posZ, channel, and type - no info in them
    dropVars = c("scanMac", "posZ", "channel", "type")
    offline = offline[ , !( names(offline) %in% dropVars ) ]
    
    # drop more unwanted access points
    offline = offline[ offline$mac %in% subMacs, ]
    
    # convert numeric values
    numVars = c("time", "posX", "posY", "orientation", "signal")
    offline[ numVars ] = lapply(offline[ numVars ], as.numeric)
    
    # convert time to POSIX
    offline$rawTime = offline$time
    offline$time = offline$time/1000
    class(offline$time) = c("POSIXt", "POSIXct")
    
    # round orientations to nearest 45
    offline$angle = roundOrientation(offline$orientation)
    
    return(offline)
  }
macs = unique(offlineSummary$mac)
online = readData("../data/online_final_trace.txt", subMacs = macs)

dist<-function(point1, point2){
  sum=0
  for (i in 1:length(point1)){
    sum=sum+(point2[i]-point1[i])^2
  }
  return(sqrt(sum))
}

knn<-function(k, data, loc, target_loc){
  temp<-data  # make a copy of the array,  not sure if needed, but being safe
  temp$dist=apply(X=temp[,-target_loc],1,FUN=dist, loc)
  temp<-temp[order(temp$dist),]
  nn=temp[1:k,target_loc]
  vote=sort(table(nn),decreasing=TRUE)[1]
  val<-names(vote)
  return(names(vote)) 
}

pdf()

tabofflineXYA = table(offline$posXY, offline$angle)
tabofflineXYA[1:6, ]

keepVars = c("posXY", "posX","posY", "orientation", "angle")
options(warn=0, error=NULL)
byLoc = with(offline,
             by(offline, list(posXY),
                function(x) {
                  ans = x[1, keepVars]
                  avgSS = tapply(x$signal, x$mac, mean)
                  y = matrix(avgSS, nrow = 1, ncol = 6,
                             dimnames = list(ans$posXY, names(avgSS)))
                  cbind(ans, y)
                }))
offlineSummary2 = do.call("rbind", byLoc)
