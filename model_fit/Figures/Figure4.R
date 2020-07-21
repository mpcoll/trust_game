library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
    "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
    "\\usepackage{amssymb}"))
tikz("ModBeh.tex", width = 6, height = 6, standAlone = TRUE,
    packages = c("\\usepackage{tikz}",
                 "\\usepackage[active,tightpage,psfixbb]{preview}",
                 "\\PreviewEnvironment{pgfpicture}",
                 "\\setlength\\PreviewBorder{0pt}",
                 "\\usepackage{amssymb}"))

par(mar=c(1,1,1,1)+0.1, mai = c(0.5,0.5,0.5,0.5), #oma = c(3,3,3,3),
cex = 1, cex.lab = 1, cex.axis =1, cex.main = 1)#, cex.text = 1)
layout(matrix(c(1,3, 2,4), 2, 2))

	Value2 = matrix(c( 0.2626 ,   0.2071,    0.2121,    0.1465 ,   0.1414 ,   0.1061,    0.1162,    0.1162,
0.0859,    0.0707, # red
0.2576,    0.2273,    0.1515,    0.1364,    0.1515,    0.1818,    0.1061,    0.0606,
0.0455,    0.0455, #brown
0.5000,    0.4697,    0.3333,    0.4091,    0.3106,    0.3409,    0.2424,    0.3106, # blue
0.2576,    0.2424,
0.4318,    0.4773,    0.3864,    0.3182,    0.3636 ,   0.2727 ,   0.3864,    0.2500, #turq
0.2500,    0.2045),10,4)



	Value5 = matrix(c(
c(0.0348,    0.0290,    0.0350,    0.0306,    0.0325,    0.0270,    0.0285,    0.0359,
0.0309 ,   0.0252),
 c( 0.0567,    0.0684,    0.0614,    0.0493 ,   0.0474 ,   0.0614 ,   0.0407,    0.0339,
0.0235   , 0.0325),
c( 0.0326,    0.0357,    0.0457,    0.0563,    0.0511,    0.0520,    0.0480,    0.0487,
0.0504,    0.0468),
 c( 0.0487,    0.0920,    0.0781,    0.0593,    0.0704,    0.0627,    0.0915,    0.0674,
0.0584,    0.0455)),10,4)

	matplot(array(c((1:10),(1:10),(1:10), (1:10) ),10,4), 
	main= "A)  Average Exchanges with perilous Trustees   ",  
	cex.main = 0.8, 
	array(c(Value2[,1], Value2[,2], Value2[,3], Value2[,4]),c(10,4)),
	axes = TRUE,
	xlab = ' ', 
	ylab = ' ', 
	xaxt = 'n',
	xlim = c(1,10), xaxp = c(1,10,9),
	 ylim=c(0,1),yaxp=c(0,1,4),
	frame.plot = FALSE, 
	type =c("l", "l", "l", "l"), lwd =2, lty = c(1,1,1,1), 
	col = c("red","coral","blue", "lightblue"), cex.axis = 1 , cex.lab = 1);
	arrows((1:10), Value2[,1], (1:10), Value2[,1]+Value5[,1], 0.01, 90);
	arrows((1:10), Value2[,1], (1:10), Value2[,1]-Value5[,1], 0.01, 90);
	arrows((1:10), Value2[,2], (1:10), Value2[,2]+Value5[,2], 0.01, 90);
	arrows((1:10), Value2[,2], (1:10), Value2[,2]-Value5[,2], 0.01, 90);
	arrows((1:10), Value2[,3], (1:10), Value2[,3]+Value5[,3], 0.01, 90);
	arrows((1:10), Value2[,3], (1:10), Value2[,3]-Value5[,3], 0.01, 90);
	arrows((1:10), Value2[,4], (1:10), Value2[,4]+Value5[,4], 0.01, 90);
	arrows((1:10), Value2[,4], (1:10), Value2[,4]-Value5[,4], 0.01, 90);
	legend(x="topright", 2, c("BPD paired Investor", "HC paired Investor", 
"BPD Trustee perilous", "HC Trustee perilous"),
       lty = c(1,1,1,1), col = 
c("blue", "lightblue","red","coral"),
	 bty = 'n', cex = 0.7)
	mtext(text = " Fraction sent ", side = 2, cex = 0.7, line = 2.2)
	mtext(text = " At Step ", side = 1, cex = 0.7, line = 2.2)
	axis(1, at= (1:10),labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), col.axis="black")




	Value9 = matrix(c( 0.3951,    0.3580 ,   0.3704 ,   0.2963,    0.3272,    0.3580,    0.3580,    0.3210,
   0.2654,    0.2531, # red
0.4091 ,   0.3939 ,   0.2955,    0.3182 ,   0.2803 ,   0.2803,    0.3182 ,   0.3030,
 0.2197,    0.2652, #brown
0.5455,    0.5795,    0.5455 ,   0.5568 ,   0.5568 ,   0.5455,    0.4432,    0.5000, # blue
 0.4205,    0.5682,
0.5741,    0.5370,    0.4907,    0.5370,    0.4815,    0.5833,    0.5463,    0.5648, #turq
0.5556,    0.4537),10,4)



	Value7 = matrix(c(
c(0.0368,    0.0341,    0.0337,    0.0348,    0.0422,    0.0395,    0.0395,    0.0397,
  0.0410,    0.0439),
c(0.0499,    0.0373,    0.0477,    0.0451,    0.0457,    0.0416,    0.0437,    0.0486,
 0.0495 ,   0.0461)),
c( 0.0396,    0.0512,    0.0490,    0.0476,    0.0596,    0.0517,    0.0612,    0.0634,
0.0672,    0.0681),
c(0.0424,    0.0530,    0.0561,    0.0677,    0.0657,    0.0650 ,   0.0568 ,   0.0678,
0.0778,    0.0701)),10,4)

	matplot(array(c((1:10),(1:10),(1:10), (1:10) ),10,4), 
	main= "B)  Average Exchanges with unperilous Trustees   ",  
	cex.main = 0.8, 
	array(c(Value9[,1], Value9[,2], Value9[,3], Value9[,4]),c(10,4)),
	axes = TRUE,
	xlab = ' ', 
	ylab = ' ', 
	xaxt = 'n',
	#xlab = 'At Step', 
	#ylab = 'Fraction sent', 
	xlim = c(1,10), xaxp = c(1,10,9),
	 ylim=c(0,1),yaxp=c(0,1,4),
	#outer = TRUE,
	frame.plot = FALSE, 
	#legend(1,2, c("LH", "HL")),
	# legend.text = c("LH", "HL"),
	# args.legend = list(x="topright", bty = 'n', cex = 2.5),
	type =c("l", "l", "l", "l"), lwd =2, lty = c(1,1,1,1), #line = 8,
	col = c("red","coral","blue", "lightblue"), cex.axis = 1 , cex.lab = 1);
	arrows((1:10), Value9[,1], (1:10), Value9[,1]+Value7[,1], 0.01, 90);
	arrows((1:10), Value9[,1], (1:10), Value9[,1]-Value7[,1], 0.01, 90);
	arrows((1:10), Value9[,2], (1:10), Value9[,2]+Value7[,2], 0.01, 90);
	arrows((1:10), Value9[,2], (1:10), Value9[,2]-Value7[,2], 0.01, 90);
	arrows((1:10), Value9[,3], (1:10), Value9[,3]+Value7[,3], 0.01, 90);
	arrows((1:10), Value9[,3], (1:10), Value9[,3]-Value7[,3], 0.01, 90);
	arrows((1:10), Value9[,4], (1:10), Value9[,4]+Value7[,4], 0.01, 90);
	arrows((1:10), Value9[,4], (1:10), Value9[,4]-Value7[,4], 0.01, 90);
	#legend(1,2, text= c("LH", "HL"), x="topright", bty = 'n', cex = 1.8)
	legend(x="topright", 2, c("BPD paired Investor", 
"HC paired Investor", 
"BPD Trustee unperilous", "HC Trustee unperilous"),
       lty = c(1,1,1,1), col = 
c("blue", "lightblue","red","coral"),
	 bty = 'n', cex = 0.7)
	mtext(text = " Fraction sent ", side = 2, cex = 0.7, line = 2.2)
	mtext(text = " At Step ", side = 1, cex = 0.7, line = 2.2)
	axis(1, at= (1:10),labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), col.axis="black")



	Value2 = matrix(c(  0.3333,    0.3421,    0.3421,    0.3333,    0.4035,    0.3684,    0.3421, # red
0.3596,    0.3070,    0.2895,
0.3356,    0.2860,    0.2523,    0.2027,    0.1847,    0.1937,    0.2050,
0.1757,    0.1284,    0.1351, #brown
0.6447,    0.6184,    0.4737,    0.6711,    0.6447,    0.7500,    0.6184, # blue
0.6579,    0.7632,    0.7632,
0.4932,    0.4899,    0.4257,    0.4189,    0.3682,    0.3750 ,   0.3378, #turq
 0.3615,    0.2838,    0.2770),10,4)



	Value5 = matrix(c(
c(0.0536,    0.0504,    0.0524,    0.0498,    0.0515,    0.0511,    0.0527,
0.0533,    0.0510,    0.0518),
c(0.1081,    0.1072,    0.1106,    0.1098,    0.1013,    0.1091,    0.1072,
0.1118,    0.1086,    0.1043),
c(0.0480,    0.0546,    0.0597,    0.0602,    0.0617,    0.0592,    0.0617,
0.0603,    0.0592,    0.0569),
 c(0.1125,    0.1289,    0.1255,    0.1304,    0.1288,    0.1233,    0.1311,
  0.1438 ,   0.1354,    0.1354)),10,4)

	matplot(array(c((1:10),(1:10),(1:10), (1:10) ),10,4), 
	main= "C)  Average Exchanges of different Risk Aversion Populations,      
$\\omega^I \\leq 1.0$ vs $\\omega^I \\geq 1.2$",  
	cex.main = 0.8, 
	array(c(Value2[,1], Value2[,2], Value2[,3], Value2[,4]),c(10,4)),
	axes = TRUE,
	xlab = ' ', 
	ylab = ' ', 
	xaxt = 'n',
	#xlab = 'At Step', 
	#ylab = 'Fraction sent', 
	xlim = c(1,10), xaxp = c(1,10,9),
	 ylim=c(0,1),yaxp=c(0,1,4),
	frame.plot = FALSE, 
	type =c("l", "l", "l", "l"), lwd =2, lty = c(1,1,1,1), #line = 8,
	col = c("red","coral","blue", "lightblue"), cex.axis = 1 , cex.lab = 1);
	arrows((1:10), Value2[,1], (1:10), Value2[,1]+Value5[,1], 0.01, 90);
	arrows((1:10), Value2[,1], (1:10), Value2[,1]-Value5[,1], 0.01, 90);
	arrows((1:10), Value2[,2], (1:10), Value2[,2]+Value5[,2], 0.01, 90);
	arrows((1:10), Value2[,2], (1:10), Value2[,2]-Value5[,2], 0.01, 90);
	arrows((1:10), Value2[,3], (1:10), Value2[,3]+Value5[,3], 0.01, 90);
	arrows((1:10), Value2[,3], (1:10), Value2[,3]-Value5[,3], 0.01, 90);
	arrows((1:10), Value2[,4], (1:10), Value2[,4]+Value5[,4], 0.01, 90);
	arrows((1:10), Value2[,4], (1:10), Value2[,4]-Value5[,4], 0.01, 90);
	#legend(1,2, text= c("LH", "HL"), x="topright", bty = 'n', cex = 1.8)
	legend(x="topright", 2, c("Investors, $\\omega^I \\leq 1.0$", "Investors, $\\omega^I \\geq 1.2$", 
"Trustees paired with $\\omega^I \\leq 1.0$", "Trustees paired with $\\omega^I \\geq 1.2$"),
       lty = c(1,1,1,1), col = 
c("blue", "lightblue","red","coral"),
	 bty = 'n', cex = 0.7)
	mtext(text = " Fraction sent ", side = 2, cex = 0.7, line = 2.2)
	mtext(text = " At Step ", side = 1, cex = 0.7, line = 2.2)
	axis(1, at= (1:10),labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), col.axis="black")

	Value =c(0.436363636363636,	0.210526315789474,
0.400000000000000,	0.500000000000000,
0.163636363636364,	0.289473684210526
)
	bar=barplot(matrix(Value,2,3), main = "D)  Trustee Guilt $(\\alpha^T)$ Distribution  ",
	cex.main = 0.8,
	xaxt = 'n',
	 ylab = " ",
	 col = c("red", "coral"),# "turquoise2", "turquoise"), 
	 axisnames = TRUE,
	 legend.text = c("BPD Trustee", "HC Trustee"),
	 args.legend = list(x="topright", bty = 'n', cex = 1),
	beside = TRUE, xaxp = c(0,3,3), ylim = c(0,1), yaxp = c(0,1,4));
	axis(1, at= (-1+3*(1:3)), labels=c("0","0.4", "1.0"), col.axis="black", lwd = 0)
	mtext(text = " Relative Empirical Frequency ", side = 2, cex = 0.7, line = 2.2)
	mtext(text = " Parameter Values ", side = 1, cex = 0.7, line = 2.2)
	text(bar[c(1 )]+0.5, 
	Value[bar[c( 2 )]]+0.3, "*", cex = 2);
	segments(bar[c(1 )], Value[bar[c( 1 )]], bar[c(1 )], Value[bar[c( 1 )]]+0.05, lwd =2.2);
	segments(bar[c(1 )], Value[bar[c( 1 )]]+0.05,  bar[c(2 )],
	 Value[bar[c( 1 )]]+0.05 , lwd=2.2);
	segments(bar[c(2 )],
	 Value[bar[c( 1 )]]+0.05 ,  bar[c(2 )],
	Value[bar[c( 2 )]], lwd=2.2);

dev.off();

tools::texi2pdf("ModBeh.tex")
system(paste(getOption("pdfviewer"), "ModBeh.pdf"))
