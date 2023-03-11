# smashburger IAQ

library(tidyverse)
library(lubridate)
library(minpack.lm)

d = read_csv('smashburger_IAQ_20230305.csv') %>%
  mutate(DATE=as.POSIXct(DATE,format="%m/%d/%Y %H:%M",tz=Sys.timezone())) %>%
  select(DATE,PM2.5,PM10,HCHO,TVOC) %>%
  pivot_longer(-DATE) %>%
  mutate(hours=as.numeric(DATE-min(DATE))/3600) %>%
  mutate(name=factor(name,levels=c('PM2.5','PM10','HCHO','TVOC'))) %>%
  arrange(name) %>%
  filter(hours<3)

ggplot(d) +
  geom_line(aes(x=hours,y=value,group=name,color=name)) +
  facet_wrap('name',scale='free_y') +
  theme_bw() +
  xlab('hours after cooking') +
  ylab('microgram per cubic meter')
ggsave('concentration_vs_time.png',width = 6,height=4)

ggplot(d) +
  geom_line(aes(x=hours,y=value,group=name,color=name)) +
  facet_wrap('name',scale='free_y') +
  theme_bw() +
  xlab('hours after cooking') +
  ylab('microgram per cubic meter') +
  scale_y_continuous(trans='log10')
ggsave('log10_concentration_vs_time.png',width = 6,height=4)




concentration = function(t=seq(0,2,by=0.05),ACH,C0,S){
  Ct = (C0- S/ACH) * exp(-t*ACH) + S/ACH
  return(Ct)
}

concentration_two_scale = function(t=seq(0,2,by=0.05),ACH=8,C0=200,S=20,lambda_slow=0.5){
  Ct = (C0- S*exp(-t*lambda_slow)/ACH) * exp(-t*ACH) + S*exp(-t*lambda_slow)/ACH
  return(Ct)
}

concentration_three_scale = function(t=seq(0,2,by=0.05),ACH=1,C0=100,S=20,lambda_slow=0.05,
                                     f1=0.5,lambda_fast=15,
                                     E=0.1,te=0.75){

  Ct = (C0- S*exp(-t*lambda_slow)/ACH) * ((1-f1)*exp(-t*ACH)+f1*exp(-t*lambda_fast)) + S*exp(-t*lambda_slow)*((1-f1)/ACH+f1/lambda_fast)

  Et = E*(C0- S*exp(-(t-te)*lambda_slow)/ACH) * ((1-f1)*exp(-(t-te)*ACH)+f1*exp(-(t-te)*lambda_fast)) + S*exp(-(t-te)*lambda_slow)*((1-f1)/ACH+f1/lambda_fast)
  Et[t<te]=0
  Ct=Ct+Et

  return(Ct)
}
plot(seq(0,10,by=0.1),log10(concentration_three_scale(t=seq(0,10,by=0.1))))


d$model=as.numeric(NA)
params=matrix(NA,nrow=4,ncol=3)
rownames(params)= unique(d$name)
colnames(params)= c('ACH','peak','source')

d$model_two_scale=as.numeric(NA)
params_two_scale=matrix(NA,nrow=4,ncol=4)
rownames(params_two_scale)= unique(d$name)
colnames(params_two_scale)= c('ACH','C0','source','lambda_slow')

d$model_three_scale=as.numeric(NA)
params_three_scale=matrix(NA,nrow=4,ncol=7)
rownames(params_three_scale)= unique(d$name)
colnames(params_three_scale)= c('ACH','C0','source','lambda_slow','f1','lambda_fast','event')


for( name in levels(d$name)){
  idx = d$name==name
  tmp=d[idx,]

  m = nlsLM(value~concentration(t=hours,ACH,C0,S),
          start=list(ACH=1,C0=max(tmp$value),S=min(tmp$value)),
          data=tmp)
  d$model[idx]=predict(m)
  params[name==unique(d$name),]=coef(m)

  m2 = nlsLM(value~concentration_two_scale(t=hours,ACH,C0,S,lambda_slow),
           start=list(ACH=1,C0=max(tmp$value),S=min(tmp$value),lambda_slow=0.1),
           data=tmp,
           lower=c(1,0,0,0),
           upper=c(Inf,Inf,Inf,1))
  d$model_two_scale[idx]=predict(m2)
  params_two_scale[name==unique(d$name),]=coef(m2)

  m3 = nlsLM(log(value)~log(concentration_three_scale(t=hours,ACH,C0,S,lambda_slow,f1,lambda_fast,E)),
             start=list(ACH=6,C0=max(tmp$value),S=min(tmp$value),lambda_slow=0.3,f1=0.5,lambda_fast=25,E=0.1),
             data=tmp,
             lower=c(0.5,0,0,0,0,15,0),
             upper=c(15,Inf,Inf,1,1,Inf,0.5),
             control=nls.lm.control(maxiter=100))
  params_three_scale[name==unique(d$name),]=coef(m3)
  d$model_three_scale[idx]=exp(predict(m3))

}

params
params_two_scale # I'm not sure why this isn't fitting reasonably, but whatever, the harder one works!
params_three_scale

ggplot(d) +
  # geom_line(aes(x=hours,y=model,group=name,color=name),size=1,linetype='dashed') +
  # geom_line(aes(x=hours,y=model_two_scale,group=name,color=name),size=1,linetype='solid') +
  geom_line(aes(x=hours,y=model_three_scale,group=name),size=1,linetype='solid') +
  geom_line(aes(x=hours,y=value,group=name,color=name)) +
  facet_wrap('name',scale='free_y') +
  scale_y_continuous(trans='log10')

# sunday's daily average PM2.5 was 4.7 mu-g/m3
# ach into the house is 1/3 per hour it looks like.

joint_model = function(t=d$hours[d$name=='PM2.5'],te=0.75,
                       lambda_slow=0.05,
                       f1=0.5,
                       lambda_fast_PM=25,
                       lambda_fast_VO=6,
                       ACH_PM=1,ACH_VO=1,
                       C0_pm2=100,S_pm2=2,E_pm=0.1,
                       C0_pm10=200,S_pm10=5,
                       C0_hcho=20,S_hcho=2,E_vo=0.5,
                       C0_tvoc=30,S_tvoc=10
                       ){

  model=c(concentration_three_scale(t,ACH=ACH_PM,lambda_slow = lambda_slow,f1=f1,lambda_fast = lambda_fast_PM,
                                    C0=C0_pm2,S=S_pm2,E=E_pm),
          concentration_three_scale(t,ACH=ACH_PM,lambda_slow = lambda_slow,f1=f1,lambda_fast = lambda_fast_PM,
                                    C0=C0_pm10,S=S_pm10,  E=E_pm),
          concentration_three_scale(t,ACH=ACH_VO,lambda_slow = lambda_slow,f1=f1,lambda_fast = lambda_fast_VO,
                                    C0=C0_hcho,S=S_hcho,  E=E_vo),
          concentration_three_scale(t,ACH=ACH_VO,lambda_slow = lambda_slow,f1=f1,lambda_fast = lambda_fast_VO,
                                    C0=C0_tvoc,S=S_tvoc,  E=E_vo)
  )
  return(model)
}


m4 = nlsLM(
           log(value)~log(joint_model(t=hours[1:(nrow(d)/4)],
                                      ACH_PM=ACH_PM, ACH_VO = ACH_VO,
                                      lambda_slow=lambda_slow,f1=f1,
                                      lambda_fast_PM=lambda_fast_PM,
                                      lambda_fast_VO=lambda_fast_VO,
                                      C0_pm2=C0_pm2,S_pm2=S_pm2,E_pm=E_pm,
                                      C0_pm10=C0_pm10,S_pm10=S_pm10,
                                      C0_hcho=C0_hcho,S_hcho=S_hcho,E_vo=E_vo,
                                      C0_tvoc=C0_tvoc,S_tvoc=S_tvoc,
                                      )),
           start=list(ACH_PM=6,ACH_VO=1,
                      lambda_slow=0.3,f1=0.5,lambda_fast_PM=20,
                      lambda_fast_VO=6,
                      C0_pm2=100,S_pm2=2,E_pm=0.1,
                      C0_pm10=200,S_pm10=5,
                      C0_hcho=1,S_hcho=0.1,E_vo=0.3,
                      C0_tvoc=1,S_tvoc=0.1
                      ),
           data=d,
           lower=c(0.5,0.5, 0,0, 15,1, rep(0,10)),
           upper=c(15,15,0.5,1,Inf,Inf,rep(Inf,10)),
           control=nls.lm.control(maxiter=100))
coef(m4)
d$model_joint=exp(predict(m4))

ggplot(d) +
  geom_line(aes(x=hours,y=model_joint,group=name),size=1,linetype='solid') +
  geom_line(aes(x=hours,y=value,group=name,color=name)) +
  facet_wrap('name',scale='free_y') +
  scale_y_continuous(trans='log10')
ggsave('log10_concentration_vs_time_joint_model.png',width = 6,height=4)

ggplot(d) +
  geom_line(aes(x=hours,y=model_joint,group=name),size=1,linetype='solid') +
  geom_line(aes(x=hours,y=value,group=name,color=name)) +
  facet_wrap('name',scale='free_y')
ggsave('concentration_vs_time_joint_model.png',width = 6,height=4)
