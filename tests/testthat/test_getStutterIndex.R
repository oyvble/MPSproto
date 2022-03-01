#we test that the getStutterIndex returns correct values
expect = function(x,y) { #helpfunction with specified tolerance
  expect_equal(as.character(x),as.character(y))
}

#library(MPSproto)

#Marker CSF1PO:
seqs = c(
"[ATCT]9",
"[ATCT]10",
"[ATCT]11",
"[ATCT]12")
type = "BW1"
x = getStutterIndex(seqs,type)
expect(x$SI,c(2,3,4,0))
expect(x$BLMM,c(10,11,12,0))

x = getStutterIndex(seqs,"BW2") #not possible
expect(x,NULL) #should return NULL

#Marker D13S317:
seqs=c(
"[TATC]11 [AATC]2 [ATCT]3 TTCT GTCT GTC",
"[TATC]12 AATC [ATCT]3 TTCT GTCT GTC",
"[TATC]12 [AATC]2 [ATCT]3 TTCT GTCT GTC",
"[TATC]13 [AATC]2 [ATCT]3 TTCT GTCT GTC")

type = "BW1"
x = getStutterIndex(seqs,type)
expect(x$SI,c(3,0,4,0))
expect(x$BLMM,c(12,0,13,0))


type = "FW1"
x = getStutterIndex(seqs,type)
expect(x$SI,c(0,0,1,3))
expect(x$BLMM,c(0,0,11,12))

#Marker D12S391:
seqs=c(
"[AGAT]8 [AGAC]7 AGAT",
"[AGAT]9 [AGAC]7 AGAT",
"[AGAT]10 [AGAC]6 AGAT",
"[AGAT]8 [AGAC]9 AGAT",
"[AGAT]10 [AGAC]7 AGAT",
"[AGAT]11 [AGAC]6 AGAT",
"[AGAT]9 [AGAC]9 AGAT",
"[AGAT]10 [AGAC]8 AGAT",
"[AGAT]10 [AGAC]5 AGAA [AGAC]3 AGAT",
"[AGAT]10 [AGAC]9 AGAT",
"[AGAT]11 [AGAC]8 AGAT",
"[AGAT]14 [AGAC]9 AGAT")

type = "BW1"
x = getStutterIndex(seqs,type)
expect(x$SI,c(2,5,6,7,0,0,10,11,0,0,0,0))
expect(x$BLMM,c(9,10,11,9,0,0,10,11,0,0,0,0))

type = "BW2" #NOTICE THE  [AGAT]9 [AGAC]9 motif
x = getStutterIndex(seqs,type)
expect(x$SI,c(0,0,5,7,8,0,0,10,0,0,0,0))
expect(x$BLMM,c(0,0,7,9,8,0,0,9,0,0,0,0))

type = "DBW1"
x = getStutterIndex(seqs,type)
expect(x$SI,c("4/5","7",0,10,0,0,0,0,0,0,0,0))
expect(x$BLMM,c("9/10",9,0,10,0,0,0,0,0,0,0,0))

type = "FWBW"
x = getStutterIndex(seqs,type)
expect(x$SI,c(0,0,2,0,0,5,0,7,0,0,10,0))
expect(x$BLMM,c(0,0,"9-7",0,0,"10-7",0,"9-9",0,0,"10-9",0))


#Marker DS21:
seqs=c(
"[TCTA]5 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11",
"[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]10",
"[TCTA]4 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11",
"[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11",
"[TCTA]4 [TCTG]7 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11",
"[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]12",
"[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]12 TA TCTA")


type = "BW1"
x = getStutterIndex(seqs,type)
expect(x$SI,c(0,4,0,6,0,0,0))
expect(x$BLMM,c(0,11,0,12,0,0,0))

type = "BW2"
x = getStutterIndex(seqs,type)
expect(x$SI,c(4,0,"4/5",0,0,0,0))
expect(x$BLMM,c("6",0,"5/7",0,0,0,0))

type = "BWFW"
x = getStutterIndex(seqs,type)
expect(x$SI,c("2",0,"1/2",0,4,0,0))
expect(x$BLMM,c("6-10",0,"5-5/5-10",0,"5-6",0,0))



