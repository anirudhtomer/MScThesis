PROC IMPORT OUT= WORK.blooddonor 
            DATAFILE= "D:\Dropbox\MSc Stats\Thesis\MScThesis\src\blood_donor\dataset.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2;
RUN;

data blooddonor;
set blooddonor;
intercep = 0.1;
run;

proc mixed data=blooddonor method=ML;
class Donate(ref="FALSE") Season(ref="Cold");
/*model Hb = intercep Age|Donate donationLast2Years|TSPD TSPD|Donate Donate|Season/solution noint;*/
model Hb = intercep Age donationLast2Years TSPD Season Donate 
donateLast2TSPD
donateLast2Donate
donateLast2Square
TSPD*TSPD
TSPD*Season

/solution noint;
random intercep donationLast2Years/type=un subject=Id vcorr=1 v=1 g gcorr;
run;


proc mixed data=blooddonor method=REML;
class Donate(ref="FALSE") Season(ref="Cold");
/*model Hb = intercep Age|Donate donationLast2Years|TSPD TSPD|Donate Donate|Season/solution noint;*/
model Hb = intercep Age donationLast2Years TSPD Season Donate 

Donate * donationLast2Years
TSPD * donationLast2Years
donationLast2Years * donationLast2Years

/solution noint;
random intercep donationLast2Years/type=un subject=Id vcorr=1 v=1 g gcorr;
run;
