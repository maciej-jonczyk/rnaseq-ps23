colorBar<-function (col, horizontal = FALSE, ...) 
{
    require(marray) || stop("Library marray is required")
    if (missing(col)) {
        col <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
            "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
            "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
            "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
            "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
            "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
            "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
            "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
            "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
            "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
            "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", 
            "#FF0060", "#FF0048", "#FF0030", "#FF0018")
    }
else {
    if (col == "fancy") {
        fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), rep(0, 
            length(c(255:150))))
        fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
        fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
            c(255:150))
        col <- rgb(b = fancy.blue/255, g = fancy.green/255, r = fancy.red/255)
    }
}
    par(mar = c(5, 2, 4, 3) + 0.1)
    maColorBar(seq(0, 1, 0.01), col = col, horizontal = FALSE, 
        k = 11, ...)
}
