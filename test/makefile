SHELL = /bin/sh

EXEC   = ../bin/rexout
DIFF   = diff -b 
MV     = mv
NCDUMP = ncdump

REX_DIR = ./rexfiles
PTS_DIR = ./pts
ANS_DIR = ./answers
BNT_DIR = ../bnt
BSN_DIR = ../sloshbsn

all: style0 style1 style2

style0: $(EXEC)
	@echo "style0"
	@$(EXEC) -rex ${REX_DIR}/camilbix.rex -pnt ${PTS_DIR}/bix.pts -out camil.bix    -style 0
	@$(EXEC) -rex ${REX_DIR}/andrwlft.rex -pnt ${PTS_DIR}/lft.pts -out andrw.lft    -style 0
	@$(EXEC) -rex ${REX_DIR}/audrybpt.rex -pnt ${PTS_DIR}/bpt.pts -out audry.bpt    -style 0
	@$(EXEC) -rex ${REX_DIR}/frederic.rex -pnt ${PTS_DIR}/emo.pts -out frederic.emo -style 0
	@$(DIFF) camil.bix    ${ANS_DIR}/style0
	@$(DIFF) andrw.lft    ${ANS_DIR}/style0
	@$(DIFF) audry.bpt    ${ANS_DIR}/style0
	@$(DIFF) frederic.emo ${ANS_DIR}/style0

style1: $(EXEC)
	@echo "style1"
	@$(EXEC) -rex ${REX_DIR}/camilbix.rex -pnt ${PTS_DIR}/bix.pts -out camil.bix    -style 1
	@$(EXEC) -rex ${REX_DIR}/andrwlft.rex -pnt ${PTS_DIR}/lft.pts -out andrw.lft    -style 1
	@$(EXEC) -rex ${REX_DIR}/audrybpt.rex -pnt ${PTS_DIR}/bpt.pts -out audry.bpt    -style 1
	@$(EXEC) -rex ${REX_DIR}/frederic.rex -pnt ${PTS_DIR}/emo.pts -out frederic.emo -style 1
	@$(DIFF) camil.bix    ${ANS_DIR}/style1
	@$(DIFF) andrw.lft    ${ANS_DIR}/style1
	@$(DIFF) audry.bpt    ${ANS_DIR}/style1
	@$(DIFF) frederic.emo ${ANS_DIR}/style1

style2: $(EXEC)
	@echo "style2"
	@$(EXEC) -rex ${REX_DIR}/camilbix.rex -pnt ${PTS_DIR}/bix.pts -bsnDir ${BSN_DIR} -style 2 -out camil.bix
	@$(EXEC) -rex ${REX_DIR}/andrwlft.rex -pnt ${PTS_DIR}/lft.pts -bsnDir ${BSN_DIR} -style 2 -out andrw.lft
	@$(EXEC) -rex ${REX_DIR}/audrybpt.rex -pnt ${PTS_DIR}/bpt.pts -bsnDir ${BSN_DIR} -style 2 -out audry.bpt
	@$(EXEC) -rex ${REX_DIR}/frederic.rex -pnt ${PTS_DIR}/emo.pts -bsnDir ${BSN_DIR} -style 2 -out frederic.emo
	@$(DIFF) camil.bix    ${ANS_DIR}/style2
	@$(DIFF) andrw.lft    ${ANS_DIR}/style2
	@$(DIFF) audry.bpt    ${ANS_DIR}/style2
	@$(DIFF) frederic.emo ${ANS_DIR}/style2

style3: $(EXEC)
	@echo "style3"
	@$(EXEC) -rex ${REX_DIR}/camilbix.rex -pnt ${PTS_DIR}/bix.pts -bsnDir ${BSN_DIR} -basin hbixdta -style 3 -out camil.bix
	@$(EXEC) -rex ${REX_DIR}/andrwlft.rex -pnt ${PTS_DIR}/lft.pts -bsnDir ${BSN_DIR} -basin lftdta  -style 3 -out andrw.lft
	@$(EXEC) -rex ${REX_DIR}/audrybpt.rex -pnt ${PTS_DIR}/bpt.pts -bsnDir ${BSN_DIR} -basin ebptdta -style 3 -out audry.bpt
	@$(EXEC) -rex ${REX_DIR}/frederic.rex -pnt ${PTS_DIR}/emo.pts -bsnDir ${BSN_DIR} -basin emobdta -style 3 -out frederic.emo
	@$(DIFF) camil.bix    ${ANS_DIR}/style3
	@$(DIFF) andrw.lft    ${ANS_DIR}/style3
	@$(DIFF) audry.bpt    ${ANS_DIR}/style3
	@$(DIFF) frederic.emo ${ANS_DIR}/style3

style4: $(EXEC)
	@echo "style4"
	@$(EXEC) -rex ${REX_DIR}/camilbix.rex -pnt ${PTS_DIR}/bix.pts -out camil.bix.nc    -bsnDir ${BSN_DIR}  -style 4 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/andrwlft.rex -pnt ${PTS_DIR}/lft.pts -out andrw.lft.nc    -bsnDir ${BSN_DIR}  -style 4 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/audrybpt.rex -pnt ${PTS_DIR}/bpt.pts -out audry.bpt.nc    -bsnDir ${BSN_DIR}  -style 4 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/frederic.rex -pnt ${PTS_DIR}/emo.pts -out frederic.emo.nc -bsnDir ${BSN_DIR}  -style 4 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(NCDUMP) camil.bix.nc    > camil.bix.nc.dump
	@$(NCDUMP) andrw.lft.nc    > andrw.lft.nc.dump
	@$(NCDUMP) audry.bpt.nc    > audry.bpt.nc.dump
	@$(NCDUMP) frederic.emo.nc > frederic.emo.nc.dump
	@$(DIFF) camil.bix.nc.dump    ${ANS_DIR}/style4
	@$(DIFF) andrw.lft.nc.dump    ${ANS_DIR}/style4
	@$(DIFF) audry.bpt.nc.dump    ${ANS_DIR}/style4
	@$(DIFF) frederic.emo.nc.dump ${ANS_DIR}/style4

style5: $(EXEC)
	@echo "style5"
	@$(EXEC) -rex ${REX_DIR}/camilbix.rex -out camil.bix.nc    -basinName hbix -bsnDir ${BSN_DIR}  -style 5 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/andrwlft.rex -out andrw.lft.nc    -basinName lft  -bsnDir ${BSN_DIR}  -style 5 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/audrybpt.rex -out audry.bpt.nc    -basinName ebpt -bsnDir ${BSN_DIR}  -style 5 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/frederic.rex -out frederic.emo.nc -basinName emob -bsnDir ${BSN_DIR}  -style 5 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(NCDUMP) camil.bix.nc    > camil.bix.nc.dump
	@$(NCDUMP) andrw.lft.nc    > andrw.lft.nc.dump
	@$(NCDUMP) audry.bpt.nc    > audry.bpt.nc.dump
	@$(NCDUMP) frederic.emo.nc > frederic.emo.nc.dump
	@$(DIFF) camil.bix.nc.dump    ${ANS_DIR}/style5
	@$(DIFF) andrw.lft.nc.dump    ${ANS_DIR}/style5
	@$(DIFF) audry.bpt.nc.dump    ${ANS_DIR}/style5
	@$(DIFF) frederic.emo.nc.dump ${ANS_DIR}/style5

style6: $(EXEC)
	@echo "style6"
	@$(EXEC) -rex ${REX_DIR}/camilbix.rex -pnt ${PTS_DIR}/bix.pts -out camil.bix.imeds    -bsnDir ${BSN_DIR}  -style 6 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/andrwlft.rex -pnt ${PTS_DIR}/lft.pts -out andrw.lft.imeds    -bsnDir ${BSN_DIR}  -style 6 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/audrybpt.rex -pnt ${PTS_DIR}/bpt.pts -out audry.bpt.imeds    -bsnDir ${BSN_DIR}  -style 6 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(EXEC) -rex ${REX_DIR}/frederic.rex -pnt ${PTS_DIR}/emo.pts -out frederic.emo.imeds -bsnDir ${BSN_DIR}  -style 6 -bnt ${BNT_DIR}/sloshdsp.bnt
	@$(DIFF) camil.bix.imeds    ${ANS_DIR}/style6
	@$(DIFF) andrw.lft.imeds    ${ANS_DIR}/style6
	@$(DIFF) audry.bpt.imeds    ${ANS_DIR}/style6
	@$(DIFF) frederic.emo.imeds ${ANS_DIR}/style6

clean:
	rm -f camil.bix*
	rm -f andrw.lft*
	rm -f audry.bpt*
	rm -f frederic.emo*
	rm -f memwatch.log
