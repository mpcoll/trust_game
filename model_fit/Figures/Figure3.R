# setPDF()
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
    "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
    "\\usepackage{amssymb}"))
tikz("ThirdCandidate.tex", width = 6, height = 9, standAlone = TRUE,
    packages = c("\\usepackage{tikz}",
                 "\\usepackage[active,tightpage,psfixbb]{preview}",
                 "\\PreviewEnvironment{pgfpicture}",
                 "\\setlength\\PreviewBorder{0pt}",
                 "\\usepackage{amssymb}"))
#pdf("BrooksSplitGraph.pdf", height = 20, width = 18, 
# bg = "white")
par(mar=c(1,1,1,1)+0.1, mai = c(0.5,0.5,0.5, 0.5), #oma = c(3,3,3,3),
cex = 1, cex.lab = 1, cex.axis =1, cex.main = 1)#, cex.text = 1)
layout(matrix(c(1,1,3,4, 1,1,3,4,2,2,5,5,2,2,5,5), 4, 4))


	Value1 =c(3.13,  13.8, 8.45 , 21.1
)
	Value1 = Value1+80;
	Value4 =c(1/sqrt(55)*43.78,  1/sqrt(38)*40.37, 1/sqrt(55)*15.7 , 1/sqrt(38)*27.6
)
	bar=barplot(matrix(Value1,2,2), main = "A)   Investment Reproduction    ",
	xaxt = 'n',
	 ylab = " ",
	 cex.axis = 1, cex.lab = 1, cex.main = 1.4,
	 col = c( "blue", "blue", "lightblue", "lightblue"),#, "cadetblue1"), 
	 density= c(-1, 12, -1, 12),
	 axisnames = TRUE,
	beside = TRUE, xaxp = c(0,5,5),#yaxt = 'n',
	 ylim = c(0,150), yaxp = c(0,150,6) )#, ylim = c(0,50), yaxp = c(0,40,6) )
	mtext(text = " Mean Investment ", side = 2, cex = 0.8, line = 2.2)
	axis(1, at= c(1.2, 2.5,  4.2, 5.5), 
	labels=c("actual", "simulated","actual","simulated"), 
	col.axis="black", lwd = 0, cex.axis = 1.2)
	mtext(text = c("BPD paired Investor", "HC paired Investor"), at = c(1.8,5),side = 1, cex = 0.8, line = 2.2)
	arrows(0.5+c((1:2), (4:5)), Value1, 0.5+c((1:2), (4:5)), Value1+Value4, 0.01, 90);
	arrows(0.5+c((1:2), (4:5)), Value1, 0.5+c((1:2), (4:5)), Value1-Value4, 0.01, 90);

	Value1 =c(18.74,   25.31 , 41.29, 37.97
)
	Value1 = Value1+60;
	Value4 =c(1/sqrt(55)*81.9,   1/sqrt(55)*24.9 , 1/sqrt(38)*101.2895, 1/sqrt(38)*30.04
)
	bar=barplot(matrix(Value1,2,2), main = "B)   Repayment Reproduction    ",
	#xlab = "Group",
	xaxt = 'n',
	  ylab = " ",
	 cex.axis = 1, cex.lab = 1, cex.main = 1.4,
	 col = c( "red", "red", "coral", "coral"),#, "cadetblue1"), 
	density= c(-1, 12, -1, 12),
	 axisnames = TRUE,
	beside = TRUE, xaxp = c(0,5,5), #yaxt = 'n',
	 ylim = c(0,150), yaxp = c(0,150,6) )#, ylim = c(0,40), yaxp = c(0,30,6), ylabels = c("80", "85", "90", "95", "100", "105", "110"));
	mtext(text = " Mean Repayment ", side = 2, cex = 0.8, line = 2.2)
	axis(1, at= c(1.2, 2.5,  4.2, 5.5), labels=c("actual", "simulated","actual","simulated"), col.axis="black", lwd = 0, cex.axis = 1.2)
	mtext(text = c("BPD Trustee", "HC Trustee"), at = c(2,5),side = 1, cex = 0.8, line = 2.2)
	arrows(0.5+c((1:2), (4:5)), Value1, 0.5+c((1:2), (4:5)), Value1+Value4, 0.01, 90);
	arrows(0.5+c((1:2), (4:5)), Value1, 0.5+c((1:2), (4:5)), Value1-Value4, 0.01, 90);

	Value5 = matrix(c(
0.5000,    0.2500,    0.2500,         0,    0.2500,    0.2500,         0, # blue
0.2500,    0.2500,    0.2500,
0.4275,    0.4288,    0.4275,    0.2163,    0.1975,    0.3700,    0.1850,
0.2225,    0.3625,    0.3300 #turq
),10,2)


	Value6 = matrix(c(
c( 0,    0,    0,    0,    0,    0,    0,
0,    0,    0),
c( 0.1711,    0.1761,    0.1747,    0.2180,    0.2250,    0.1643,    0.2245,
0.2182,    0.1519,    0.1296)),10,2)


	matplot(array(c((1:10),(1:10)),10,2), array(c(Value5[,1], Value5[,2]),c(10,2)),
	main= "C)     Reproduced Investor Trajectory         ",  
	cex.main = 1.4, 	
	#axes = TRUE,
	yaxt = 't',
	xlab = ' ',#'At Step', 
	ylab = ' ', 
	xlim = c(1,10), #xaxp = c(1,10,9),
	xaxt = 'n',
	 ylim=c(0,1),yaxp=c(0,1,4),
	frame.plot = FALSE, 
	type =c( "l", "l"), lwd = c(3,1), lty = c(1,2), #line = 8,
	col = c("blue","blue"), cex.axis = 1 , cex.lab = 1);
	polygon( c((1:10),  (10:1)), c(pmin(Value5[1:10,2]+Value6[1:10,2],1.0), 
	pmax(Value5[10:1,2]-Value6[10:1,2],0.0)), 
	col = rgb(0,0,200, alpha =30, maxColorValue=255), density = 20, 
	border = 'blue', lty = 2);
	mtext(text = " Fraction sent ", side = 2, cex = 0.8, line = 2.2)
	mtext(text = " At Step ", side = 1, cex = 0.8, line = 2.2)
	axis(1, at= (1:10),labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), col.axis="black")
	legend(x="topright", 2, c("Investor", "Simulated Investor"),
 lwd =c(3,1), lty = c(1,2),
      #fill =c("blue","blue",  "red","red"), 
	col = c("blue","blue"), #density = c(-1, 20, -1, 20),
	 bty = 'n', cex = 1)

	Value5 = matrix(c(  0.5000,    0.1667,         0 ,        0,    0.1667,         0,         0,  # red
 0.6667,    0.3333,         0,
0.2592,    0.1825,    0.1275 ,        0 ,   0.1483,    0.1342,         0,
 0.2142,    0.1150,    0.0542 # coral
),10,2)


	Value6 = matrix(c(
c( 0,    0,    0,    0,    0,    0,    0,
   0,    0,    0),
 c( 0.1375,    0.1757,    0.1679,         0,    0.1866,    0.1501,         0,
    0.1908,    0.1507,    0.1029)),10,2)


	matplot(array(c((1:10),(1:10)),10,2), array(c(Value5[,1], Value5[,2]),c(10,2)),
	main= "D)     Reproduced Trustee Trajectory         ",  
	cex.main = 1.4, 	
	#axes = TRUE,
	yaxt = 't',
	xlab = ' ',#'At Step', 
	ylab = ' ', 
	xlim = c(1,10), #xaxp = c(1,10,9),
	xaxt = 'n',
	 ylim=c(0,1),yaxp=c(0,1,4),
	#outer = TRUE,
	frame.plot = FALSE, 
	#legend(1,2, c("LH", "HL")),
	# legend.text = c("LH", "HL"),
	# args.legend = list(x="topright", bty = 'n', cex = 2.5),
	type =c("l", "l"), lwd =c(3,1), lty = c(1,2), #line = 8,
	col = c("red", "red"), cex.axis = 1 , cex.lab = 1);
	polygon( c((1:10),  (10:1)), c(pmin(Value5[1:10,2]+Value6[1:10,2],1.0), 
	pmax(Value5[10:1,2]-Value6[10:1,2],0.0)), 
	col = rgb(200,0,0, alpha =30, maxColorValue=255), density = 20, 
	border = 'red', lty = 2);
	mtext(text = " Fraction sent ", side = 2, cex = 0.8, line = 2.2)
	mtext(text = " At Step ", side = 1, cex = 0.8, line = 2.2)
	axis(1, at= (1:10),labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), col.axis="black")
	legend(x="topright", 2, c("Trustee", "Simulated Trustee"), 
lwd =c(3,1), lty = c(1,2),
      #fill =c("blue","blue",  "red","red"), 
	col = c("red", "red"), #density = c(-1, 20, -1, 20),
	 bty = 'n', cex = 1)

dev.off();

tools::texi2pdf("ThirdCandidate.tex")
system(paste(getOption("pdfviewer"), "ThirdCandidate.pdf"))
