library(MASS)
m=read.delim("s5g.tsv"); m$fsm=as.integer(m$clade_x=="FSM"); mag=m[m$kind=="MAG"&m$clade_x%in%c("FSM","RGM"),]
for(f in c("pg~fsm","pg~fsm+cpm","pg~fsm+cpm+completeness_used_pct","pg~fsm+cpm+checkm2_this_study_completeness_pct","pg~fsm+cpm+completeness_pct")){
 try({r=glm.nb(as.formula(f),data=mag);cat(f,coef(r)["fsm"],summary(r)$coefficients["fsm",4],"\n")})}
for(f in c("pg~fsm+offset(log(mb))","pg~fsm+cpm+offset(log(mb))","pg~fsm+cpm+completeness_used_pct+offset(log(mb))","pg~fsm+cpm+checkm2_this_study_completeness_pct+offset(log(mb))")){
 try({r=glm.nb(as.formula(f),data=mag);cat(f,coef(r)["fsm"],summary(r)$coefficients["fsm",4],"\n")})}
p=m[m$clade_x%in%c("FSM","RGM"),]; 
for(f in c("pg~fsm+offset(log(mb))","pg~fsm+cpm+offset(log(mb))")){r=glm.nb(as.formula(f),data=p);cat("325",nrow(p),f,coef(r)["fsm"],exp(coef(r)["fsm"]),summary(r)$coefficients["fsm",4],"\n")}
