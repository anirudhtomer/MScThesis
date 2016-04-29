PROC IMPORT OUT= WORK.blooddonor 
            DATAFILE= "D:\Dropbox\MSc Stats\Thesis\MScThesis\latex\code\blood_donor\dataset.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2;
RUN;

data blooddonor;
set blooddonor;
intercep = 0.1;
/*donationLast2Years_10 = donationLast2Years/100;*/
run;

proc mixed data=blooddonor method=ML;
class Donate(ref="FALSE") Season(ref="Cold");
/*model Hb = intercep Age|Donate donationLast2Years|TSPD TSPD|Donate Donate|Season/solution noint;*/
model Hb = intercep Age donationLast2Years TSPD Season Donate 

Age * Season
donationLast2Years * donationLast2Years
donationLast2Years * TSPD
donationLast2Years * Season
donationLast2Years * Donate
TSPD * Season
TSPD * Donate


/solution noint;
random intercep donationLast2Years/type=un subject=Id vcorr=1 v=1 g gcorr;
run;
