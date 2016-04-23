PROC IMPORT OUT= WORK.blooddonor 
            DATAFILE= "D:\Dropbox\MSc Stats\Thesis\MScThesis\latex\code\blood_donor\dataset.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2;
RUN;

data blooddonor;
set blooddonor;
intercep = 0.1;
donationLastTwoYears_10 = donationLastTwoYears/100;
run;

proc mixed data=blooddonor method=REML;
class Donate(ref="FALSE") Season(ref="Cold");
model Hb = intercep Donate Season donationLastTwoYears_10 Age TSPD/solution noint;
random intercep donationLastTwoYears_10/type=un subject=Id vcorr=1 v=1 g gcorr;
run;
