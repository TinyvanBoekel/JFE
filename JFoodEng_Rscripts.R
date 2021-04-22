
# R script for the manuscript "Kinetics of heat-induced chnages in foods: a workflow proposal" in J. Food Eng. by M.A.J.S.van Boekel

# libraries used:
library(papaja)
library(rmarkdown)
library(bookdown)
library(ggplot2)
library(dplyr)
library(brms)
library(tidyverse)
library(rstan)
library(broom)
library(GGally)
library(tidybayes)
library(here)
library(patchwork)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())


knitr::opts_chunk$set(echo = FALSE, 
                      fig.align = "center",
                      fig.width = 6, 
                      fig.height = 4, 
                      dev = "png",
                    warning = FALSE,
                    message = FALSE,
                    results= "asis",
                    comment = NA)

theme_set(theme_bw())

# import of data:
carnitin <- read.csv(here("data", "carnitine_data.csv"), header=TRUE, sep=",")

# rescale time from minutes to hours
carnitin <- carnitin %>%  mutate(time=time/60)
#preparing the data at each temp:
carn353 <- subset(carnitin, temp ==353, select=time:conc)
carn358 <- subset(carnitin, temp ==358, select=time:conc)
carn363 <- subset(carnitin, temp ==363, select=time:conc)
carn368 <- subset(carnitin, temp ==368, select=time:conc)
carn373 <- subset(carnitin, temp ==373, select=time:conc)
carn383 <- subset(carnitin, temp ==383, select=time:conc)
carn393 <- subset(carnitin, temp ==393, select=time:conc)
carn403 <- subset(carnitin, temp ==403, select=time:conc)

# plotting the data:
carnitin %>% mutate(temp=str_c("T=",temp)) %>% ggplot(aes(x=time, y=conc, group=temp))+
  geom_point(aes(y=conc))+
  facet_wrap(~temp, ncol=4, scales = "free")+
  labs(x = "time (h)", y = "carnitin (mg/L)")+
  theme(strip.background = element_rect(color="black", fill="lightblue", size=1.5, linetype="solid")   )


# nonlinear individual regressions with brms:

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), 
           c0~1, 
           k~1, 
           nt~1,
           nl=TRUE)

nlprior<-c(prior(normal(66,10), nlpar = "c0"),
              prior(normal(1,5), nlpar="k", lb=0),
              prior(normal(1,1), nlpar="nt"),
           prior(cauchy(0,5), class="sigma")
           )           
 
model353 <- brm(formula=nlform, data=carn353, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits", "fit353"))

model358 <- brm(formula=nlform, data=carn358, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit358"))

model363 <- brm(formula=nlform, data=carn363, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit363"))

model368 <- brm(formula=nlform, data=carn368, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit368"))

model373 <- brm(formula=nlform, data=carn373, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit373"))

model383 <- brm(formula=nlform, data=carn383, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15,file=here("fits","fit383"))

model393 <- brm(formula=nlform, data=carn393, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.999), seed=15, file=here("fits","fit393"))

model403 <- brm(formula=nlform, data=carn403, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.999), seed=15, file=here("fits","fit403"))

# collecting individual parameter estimates:
post1 <- posterior_samples(model353)
post2 <- posterior_samples(model358)
post3 <- posterior_samples(model363)
post4 <- posterior_samples(model368)
post5 <- posterior_samples(model373)
post6 <- posterior_samples(model383)
post7 <- posterior_samples(model393)
post8 <- posterior_samples(model403)

p <-
  bind_rows(
  post1,
  post2,
  post3,
  post4,
  post5,
  post6,
  post7,
  post8
)
iter <- 8000

# plotting parameter estimates from individual regressions:
p <- 
  p %>% 
  mutate(temperature = rep(c("353 K","358 K","363 K","368 K","373 K","383 K","393 K","403 K"),each = iter)) %>% mutate(log_kr=log(b_k_Intercept))

plot_nt <- p %>% 
  ggplot(aes(x = b_nt_Intercept, y = temperature)) +
  geom_halfeyeh(fill = "lightblue", point_interval = mean_qi, .width = .95)+ 
  labs(x=expression(n[t]))

plot_c0 <- p %>% 
  ggplot(aes(x = b_c0_Intercept, y = temperature)) +
  geom_halfeyeh(fill = "lightblue", 
                point_interval = mean_qi, .width = .95) +
  labs(x=expression(c[0]), y="")

plot_kr <- p %>% 
  ggplot(aes(x = log_kr, y = temperature)) +
  geom_halfeyeh(fill = "lightblue", point_interval = mean_qi, .width = .95) +
  labs(x=expression(ln(k[r])), y="")

plot_c0+plot_nt+plot_kr

# example correlation plot at T=373:

cor_373 <- dplyr::select(post5,b_c0_Intercept:sigma)
# change the names of the columns to be displayed in the panels
cor_373 <- setNames(cor_373, c(expression(paste("c"[0])), expression(paste("k"[r])), expression(paste("n"[t])),expression(sigma[e])))

# use ggally for a pairs plot
(cor_plot373 <-cor_373  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
           mapping=aes(fill="red"),
           upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef")),
           labeller=label_parsed)+ 
          theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")
          ))

# brms multilevel model for nt:
                   
nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*kref*exp(Ea*(1-373.5747/temp))*time)^(1/(1-nt)), 
           c0~1,
           nt~1+(1|temp),
           kref~1,
           Ea ~1, 
           nl=TRUE)

# set priors for the parameters:
nlprior<-c(prior(normal(67,5), nlpar = "c0", lb=0),
              prior(normal(0,1), nlpar="kref"),
                prior(normal(1,0.5), nlpar="nt"),
                 prior(normal(27, 1), nlpar="Ea"),
           prior(cauchy(0,25), class="sigma")
          )

model_all_multi<-brm(formula=nlform, data=carnitin, family = gaussian(), prior = nlprior, warmup=2000, iter=4000, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits","model_all_multi"))

# transferring the posterior to a dataframe:
model_all_multi_post <- posterior_samples(model_all_multi)

#print(model_all_multi, digits=4)
#pairs(model_all_multi)

# correlation plot using GGally:

cor_alldata <- dplyr::select(model_all_multi_post,b_c0_Intercept:sigma)
# change the names of the columns to be displayed in the panels
cor_alldata <- setNames(cor_alldata, c(expression(paste("c"[0])), expression(paste("n"[t])), expression(paste("k"[ref])), expression(paste("E"[a],"/RT"[ref])), expression(sigma[n[t]]), expression(sigma[e])))

(cor_alldata_plot <-cor_alldata  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
           mapping=aes(fill="red"),
           upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef.")), 
           labeller=label_parsed)+ 
          theme(strip.text.x = element_text(size = 10, color = "red"),
                    strip.text.y = element_text(size = 10, color = "red")) +
                       theme(axis.text.x = element_text(angle=70, hjust=1))
          )

# back transform Ea and preparing for table using apa_table from papaja:
p1 <- summary(model_all_multi)
summary_p1 <- rbind(data.frame(p1$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(p1$random$temp) %>%  dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(p1$spec_pars)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

summary_p1[4,] <- summary_p1[4,]*8.314*373.5747/1000

rownames(summary_p1) <- c("$c_0 \\text { (g/dm}^3)$", "$n_t \\text { (-)}$", "$k_{ref} \\text { ((dm}^3 \\text{mol}^{-1})^{n_t-1} \\text {h}^{-1})$","$E_a \\text{ (kJ/mol)}$","$\\sigma_{n_t}$", "$\\sigma_e$")
colnames(summary_p1) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
    summary_p1,
    placement = "H",
    align = c("c", "c", "c", "c"),
    caption = "(ref:all-data-tab)",
    note = NULL,
    small = TRUE,
    format.args=list(
      digits = c(2,2,2,2)
    ),
    escape = FALSE
    )

# plot of random parameter nt:

mcmc_plot(model_all_multi, pars=c("^r_"), point_est="mean", prob_outer=0.95, prob=0.5) + ggplot2::scale_y_discrete(labels=c(expression(paste("n"[t], " , T=353")),expression(paste("n"[t], " , T=358")),expression(paste("n"[t], " , T=363")),expression(paste("n"[t], " , T=368")),expression(paste("n"[t], " , T=373")),expression(paste("n"[t], " , T=383")),expression(paste("n"[t], " , T=393")) ,expression(paste("n"[t], " , T=403"))))+
  labs(x=expression(paste("deviations of n"[t]," at the group level from the population level n"[t], "=1.14")))

# combine the data with predictions using the random effects with re_formula = NULL
newvary1 <- expand.grid(time=seq(from = 0, to =10, by=0.5),temp=c(353, 358))
newvary2 <- expand.grid(time=seq(from = 0, to =6, by=0.25),temp=c(363, 368, 373))
newvary3 <- expand.grid(time=seq(from = 0, to =3, by=0.1),temp=c(383, 393))
newvary4 <- expand.grid(time=seq(from = 0, to =1, by=0.05),temp=c(403))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4)

carn_ind <- cbind(newvary, predict(model_all_multi, newdata=newvary, re_formula = NA)[,-2])  # fit for the population level
carn_ind2 <- cbind(newvary, predict(model_all_multi, newdata=newvary, re_formula = NULL[,-2])) #fit for the group level
names(carn_ind) <- c("time", "temp", "conc", "lower", "upper")

carn_ind$temp=as.factor(carn_ind$temp)

(ind_fit <-  ggplot(carnitin, aes(x=time, y=conc))+
  geom_point(colour = "#2c3e50",
    fill = "#2c3e50")+
  facet_wrap(~temp, ncol=4, scales="free_x")+
  geom_line(data = carn_ind, aes(y = conc), size = 1, colour="blue") +
    geom_line(data = carn_ind, aes(y = lower), lty = 2) +
    geom_line(data = carn_ind, aes(y = upper), lty = 2)) +
  geom_line(data = carn_ind2, aes(y = Estimate), size = 1, colour="red") + #fit for the population level
    labs(x="time (h)", y="carnitin (mg/L)")

#The following brms code is to analyze completely pooled data, so it is not multilevel but complete pooling

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*kref*exp(Ea*(1-373.5747/temp))*time)^(1/(1-nt)), 
           c0~1,
           nt~1,
           kref~1,
           Ea ~1, 
           nl=TRUE)

# set priors for the parameters:
nlprior<-c(prior(normal(67,5), nlpar = "c0", lb=0),
              prior(normal(0,1), nlpar="kref"),
                prior(normal(1,0.5), nlpar="nt"),
                 prior(normal(27, 1), nlpar="Ea"),
           prior(cauchy(0,25), class="sigma")
          )

model_all_nth2<-brm(formula=nlform, data=carnitin, family = gaussian(), prior = nlprior, warmup=2000, iter=4000, control = list(adapt_delta = 0.99, max_treedepth=15), file=here("fits","model_all_nth2"))
#print(model_all_nth2, digits=4)

#Tb = 373.0, code for first-order model equation:
nlform<-bf(conc ~ c0*exp(-kref*exp(Ea*(1-373.0/temp))*time), 
           c0+kref+Ea ~1, 
           nl=TRUE)

# set priors for the parameters:
nlprior<-c(prior(normal(66,1), nlpar = "c0"),
              prior(normal(0,2), nlpar="kref"),
               prior(normal(27, 1), nlpar="Ea"),
                  prior(cauchy(0,25), class="sigma")
           )

model_all_fo<-brm(formula=nlform, data=carnitin, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, control = list(adapt_delta = 0.99), seed=15, file=here("fits","model_all_fo"))
model_all_fo_post <- posterior_samples(model_all_fo)
#print(model_all_fo, digits=4)

# preparing for papaja table:
pp1 <- summary(model_all_fo)
summary_pp1 <- rbind(data.frame(pp1$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(pp1$spec_pars)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

# back transform Ea:
summary_pp1[3,] <- summary_pp1[3,]*8.314*373.5747/1000

rownames(summary_pp1) <- c("$c_0 \\text { (g/dm}^3)$","$k_{ref} \\text { (h}^{-1})$","$E_a \\text{ (kJ/mol)}$", "$\\sigma_e$")
colnames(summary_pp1) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
    summary_pp1,
    placement = "H",
    align = c("c", "c", "c", "c"),
    caption = "(ref:alldata-fo-summary)",
    note = NULL,
    small = TRUE,
    format.args=list(
      digits = c(2,2,2,2)
    ),
    escape = FALSE
    )

# preparing for loo-cv and reporting with papaja in an apa_table:
model_all_multi <- add_criterion(model_all_multi, c("loo"), file="model_all_multi")
model_all_fo <- add_criterion(model_all_fo, c("loo"), file="model_all_fo")
model_all_nth2 <- add_criterion(model_all_nth2, c("loo"), file="model_all_nth2")

loo_nt <- loo(model_all_multi)
loo_fo <- loo(model_all_fo)
loo_nth2 <- loo(model_all_nth2)
loo_result <- loo_compare(loo_nt,loo_nth2,loo_fo)

elpd_diff <- c(loo_result[1,1], loo_result[2,1], loo_result[3,1])
se_diff <- c(loo_result[1,2], loo_result[2,2], loo_result[3,2])
loo_df <- data.frame(elpd_diff, se_diff)
colnames(loo_df) <- c("elpd-difference","SE" )
rownames(loo_df) <- c("multilevel n-th order", "single-level n-th order", "single level first-order")

apa_table(
  loo_df,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:loo-comp)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)
# posterior predictive checks:
ppcA <- pp_check(model_all_multi, nsamples = 100) + labs(x="carnitin, mg/L",subtitle = "A")
ppcC <- pp_check(model_all_fo, nsamples = 100) + labs(x="carnitin, mg/L",subtitle = "C")
ppcB <- pp_check(model_all_nth2, nsamples=100) +labs(x="carnitin, mg/L", subtitle="B")
ppcA + ppcB + ppcC

# Codes for the classical two-step approach using least-squares linear regression:
# 80 C/353 K:

ln_carn353 <- transform(carn353,conc=log(conc))
# lm result (frequentist approach):
m353 <- lm(conc~1+time, data=ln_carn353)
summary(lm353)

# bayesian approach with brms:
carn353fit <- 
  brm(data = ln_carn353, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(0.4, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","carn353fit"))
# 85 C/358 K:

ln_carn358 <- transform(carn358,conc=log(conc))
# lm result:
lm358 <- lm(conc~1+time, data=ln_carn358)
summary(lm358)
#linearized first order model

carn358fit <- 
  brm(data = ln_carn358, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(0.4, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","carn358fit"))

# 90 C/363 K:

ln_carn363 <- transform(carn363,conc=log(conc))

# lm result:
lm363 <- lm(conc~1+time, data=ln_carn363)
summary(lm363)

carn363fit <- 
  brm(data = ln_carn363, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(0.6, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","carn363fit"))

# 95 C/368 K:

ln_carn368 <- transform(carn368,conc=log(conc))
carn368fit <- 
  brm(data = ln_carn368, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(0.6, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","carn368fit"))

# 100 C/373 K:

ln_carn373 <- transform(carn373,conc=log(conc))
carn373fit <- 
  brm(data = ln_carn373, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(0.8, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","carn373fit"))

# 110/383 K:

ln_carn383 <- transform(carn383,conc=log(conc))
carn383fit <- 
  brm(data = ln_carn383, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(0.8, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95),file=here("fits","carn383fit"))

# 120/393 K:

ln_carn393 <- transform(carn393,conc=log(conc))
carn393fit <- 
  brm(data = ln_carn393, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(1, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","carn393fit"))

# 130/403 K:

ln_carn403 <- transform(carn403,conc=log(conc))
carn403fit <- 
  brm(data = ln_carn403, family = gaussian,
      formula = conc ~ 1 + time,
      prior = c(set_prior("normal(4.6, 50)", class = "Intercept"),
                set_prior("normal(1, 1)", class = "b"),
                set_prior("cauchy(0,25)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","carn403fit"))

# preparing for Arrhenius fit:
Tinv <- c(1/(8.314*353),1/(8.314*358),1/(8.314*363),1/(8.314*368),1/(8.314*373),1/(8.314*383),1/(8.314*393),1/(8.314*403))
ln_k <- c(log(-fixef(carn353fit)[2,1]),log(-fixef(carn358fit)[2,1]),log(-fixef(carn363fit)[2,1]), log(-fixef(carn368fit)[2,1]),log(-fixef(carn373fit)[2,1]), log(-fixef(carn383fit)[2,1]),log(-fixef(carn393fit)[2,1]),log(-fixef(carn403fit)[2,1]))
df_linArrh <- data.frame(Tinv,ln_k)

linArrhfit <- 
  brm(data = df_linArrh, family = gaussian,
      formula = ln_k ~ 1 + Tinv,
      prior = c(set_prior("normal(25, 10)", class = "Intercept"),
                set_prior("normal(-100000, 10000)", class = "b"),
                set_prior("cauchy(0,10)", class = "sigma")),
      chains = 4, iter = 4000, warmup = 2000, control = list(adapt_delta =0.95), file=here("fits","linArrhfit"))

# calculate k0 from its ln, and recalculate Ea in kJ/mol:
linArrh_post <- posterior_samples(linArrhfit) %>% mutate(b_Tinv=-b_Tinv) %>% mutate(k0=exp(b_Intercept)) %>% mutate(Ea=b_Tinv/1000)

# range to plot:
T.seq <- data.frame(Tinv = seq(from = 1/(8.314*405), to = 1/(350*8.314), by = 0.00001/8.314))

#fitted is about mu:
muSummary <-
  fitted(linArrhfit, 
         newdata = T.seq) %>%
  as_tibble() %>%
  bind_cols(T.seq)

#predict is about future individual values:
pred.Arrh <-
  predict(linArrhfit,
          newdata = T.seq) %>%
  as_tibble() %>%
  bind_cols(T.seq)

#plot of fitted and predicted values
plot_Arrh <- df_linArrh %>%
  ggplot(aes(x = Tinv, y = ln_k)) +
  geom_ribbon(data = pred.Arrh, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = muSummary, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = muSummary, aes(y = Estimate)) +
   geom_point(color = "navyblue", shape = 16, size = 1.5, alpha = 2/3)+
  labs(x="1/RT", y="ln k", subtitle="A")+
   theme(axis.text.x = element_text(angle=70, hjust=1))

resid1 = resid(linArrhfit)[, "Estimate"]
residplot1 <-  ggplot() + geom_point(data = NULL, size = 2, aes(y = resid1, 
                                                              x = df_linArrh$Tinv))+
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  labs(x = "1/RT", y = "residual ln k", subtitle = "B")+
   theme(axis.text.x = element_text(angle=70, hjust=1))

# calculating QQ plot for Arrhenius residuals:

res_linArrh <- as.data.frame(residuals(linArrhfit, summary=T))
qqplot_linArrh <-  ggplot(res_linArrh, aes(sample = Estimate)) +
  stat_qq() +
  stat_qq_line()+
  labs(subtitle = "C")+
   theme(axis.text.x = element_text(angle=70, hjust=1))

plot_Arrh + residplot1 + qqplot_linArrh

# correlation plot for Arrhenius:
linArrh_cor <- dplyr::select(linArrh_post,b_Intercept:sigma)
# change the names of the columns to be displayed in the panels
linArrh_cor <- setNames(linArrh_cor, c(expression(paste("ln k"[0])), 
      expression(paste("E"[a])), expression(sigma[e])))

# use ggally for a pairs plot
(linArrh_cor_plot <-linArrh_cor  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
           mapping=aes(fill="red"),
           upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr.coef.")), 
           labeller=label_parsed)+ 
           theme(strip.text.x = element_text(size = 10, color = "red"),
                  strip.text.y = element_text(size = 10, color = "red")) +
           theme(axis.text.x = element_text(angle=70, hjust=1))         
  )

# reporting Arrhenius results in papaja table:
fit_summary1 <- summary(linArrhfit)
fit_summary2 <- rbind(data.frame(fit_summary1$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(fit_summary1$spec_pars)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))
fit_summary2[2,1]=-fit_summary2[2,1]
fit_summary2[2,3]=-fit_summary2[2,3]
fit_summary2[2,4]=-fit_summary2[2,4]
fit_summary2[2,]=fit_summary2[2,]/1000
rownames(fit_summary2) <- c("$\\ln k_0 \\text { (h}^{-1})$","$E_a \\text{ (kJ/mol)}$", "$\\sigma_e$")
colnames(fit_summary2) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
    fit_summary2,
    placement = "H",
    align = c("c", "c", "c","c"),
    caption = "(ref:linArrh-table)",
    note = NULL,
    small = TRUE,
    format.args=list(
      digits = c(2,2,2,2)
    ),
    escape = FALSE
    )

# preparing to plot Ea and lnk0 densities:
model_all_fo_post <- model_all_fo_post %>% mutate(Ea=b_Ea_Intercept*8.314*373.5747/1000) %>% mutate(lnk0=log(b_kref_Intercept/(exp(-Ea*1000/(8.314*373.57))))) 

Ea_density <- ggplot(data=linArrh_post, aes(x=Ea))+
  geom_density(fill="red", alpha=0.4)+ 
  geom_density(data=model_all_fo_post, aes(x=Ea), fill="blue", alpha=0.3)+
  labs(x=expression(paste("E"[a], " (kJ/mol)")), subtitle = "A")  

lnk0_density <- ggplot(data=linArrh_post, aes(x=b_Intercept))+
  geom_density(fill="red", alpha=0.4)+
  geom_density(data=model_all_fo_post, aes(x=lnk0), fill="blue",alpha=0.3)+
  labs(x=expression(paste("ln k"[0])," 1/s"), subtitle = "B")

Ea_density + lnk0_density

# nonlinear bayesian regression for first-order model:
nlform <- bf(conc~ c0*exp(-kr*time), c0~1, kr~1, nl=TRUE)
nlprior <- c(prior(normal(66,10),nlpar="c0" ),
             prior(normal(1,0.5), nlpar="kr")
)
model353_fo <- brm(formula=nlform, data=carn353, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model353_fo"))

model353_fo_post <- posterior_samples(model353_fo)

model358_fo <- brm(formula=nlform, data=carn358, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model358_fo"))

model358_fo_post <- posterior_samples(model358_fo)

model363_fo <- brm(formula=nlform, data=carn363, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model363_fo"))

model363_fo_post <- posterior_samples(model363_fo)

model368_fo <- brm(formula=nlform, data=carn368, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model368_fo"))

model368_fo_post <- posterior_samples(model368_fo)

model373_fo <- brm(formula=nlform, data=carn373, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model373_fo"))

model373_fo_post <- posterior_samples(model373_fo)

model383_fo <- brm(formula=nlform, data=carn383, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model383_fo"))

model383_fo_post <- posterior_samples(model383_fo)

nlform <- bf(conc~ c0*exp(-kr*time), c0~1, kr~1, nl=TRUE)
nlprior <- c(prior(normal(66,10),nlpar="c0" ),
             prior(normal(3,2), nlpar="kr")
)

model393_fo <- brm(formula=nlform, data=carn393, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model393_fo"))

model393_fo_post <- posterior_samples(model393_fo)

model403_fo <- brm(formula=nlform, data=carn403, family=gaussian(),
                   prior=nlprior, control = list(adapt_delta=0.99),
                   file=here("fits","model403_fo"))

model403_fo_post <- posterior_samples(model403_fo)

#compiling rate constants from linear regression:
kr_lr <- c(-fixef(carn353fit)[2,1],-fixef(carn358fit)[2,1],-fixef(carn363fit)[2,1], -fixef(carn368fit)[2,1],-fixef(carn373fit)[2,1], -fixef(carn383fit)[2,1],-fixef(carn393fit)[2,1],-fixef(carn403fit)[2,1])

#compiling rate constants from nonlinear regression:
kr_nlr <- c(fixef(model353_fo)[2,1],fixef(model358_fo)[2,1],fixef(model363_fo)[2,1], fixef(model368_fo)[2,1],fixef(model373_fo)[2,1], fixef(model383_fo)[2,1],fixef(model393_fo)[2,1],fixef(model403_fo)[2,1])

# plotting comparison
df_kr_lr <- data.frame(kr_lr, kr_nlr)
df_kr_lr %>% ggplot(aes(x=kr_lr, y=kr_nlr))+
  geom_point()+
  geom_abline(intercept = 0, lty=2)+
  labs(x="linear regression estimate", y="nonlinear regression estimate")


# prediction of k
linArrh_post <- linArrh_post %>% mutate(k343=exp(b_Intercept-Ea*1000/(343*8.314))) 
model_all_fo_post <- model_all_fo_post %>% mutate(k343=b_kref_Intercept*exp(b_Ea_Intercept*(1-373.0/343.0)))

k343_density <- ggplot(data=linArrh_post, aes(x=k343))+
  geom_density(fill="turquoise", alpha=0.4)+ 
  geom_density(data=model_all_fo_post, aes(x=k343), fill="green", alpha=0.6)+
  labs(x=expression(paste(k[343]))) 


linArrh_post <- linArrh_post %>% mutate(k413=exp(b_Intercept-Ea*1000/(413*8.314))) 
model_all_fo_post <- model_all_fo_post %>% mutate(k413=b_kref_Intercept*exp(b_Ea_Intercept*(1-373.0/413.0)))

k413_density <- ggplot(data=linArrh_post, aes(x=k413))+
  geom_density(fill="turquoise", alpha=0.4)+ 
  geom_density(data=model_all_fo_post, aes(x=k413), fill="green", alpha=0.6)+
  labs(x=expression(paste(k[413]))) 

k343_density + k413_density


# The following code is used in the Supplement

# prior predictive check with brms:

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), 
           c0~1, 
           k~1, 
           nt~1,
           nl=TRUE)

nlprior<-c(prior(normal(66,2), nlpar = "c0"),
           prior(normal(1,5), nlpar="k", lb=0),
           prior(normal(1,0.5), nlpar="nt", lb=0),
           prior(cauchy(0,5), class="sigma")
)           

model353a <- brm(formula=nlform, data=carn353, family = gaussian(), prior = nlprior, sample_prior="only", iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit353a"))

plot(conditional_effects(model353a, spaghetti=TRUE, nsamples=500))

# trace plot for result at T=353

model353_post <- posterior_samples(model353) %>% select(-lp__)
bayesplot::mcmc_trace(model353_post)

# code for individual fits:

nlform<-bf(conc ~ (c0^(1-nt)+(nt-1)*k*time)^(1/(1-nt)), 
           c0~1, 
           k~1, 
           nt~1,
           nl=TRUE)

nlprior<-c(prior(normal(66,2), nlpar = "c0"),
           prior(normal(1,5), nlpar="k", lb=0),
           prior(normal(1,0.5), nlpar="nt"),
           prior(cauchy(0,5), class="sigma")
)           

model353 <- brm(formula=nlform, data=carn353, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit353"))

model358 <- brm(formula=nlform, data=carn358, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit358"))

model363 <- brm(formula=nlform, data=carn363, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file="fit363")

model368 <- brm(formula=nlform, data=carn368, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit368"))

model373 <- brm(formula=nlform, data=carn373, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit373"))

model383 <- brm(formula=nlform, data=carn383, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.99), seed=15, file=here("fits","fit383"))

model393 <- brm(formula=nlform, data=carn393, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.999), seed=15, file=here("fits","fit393"))

model403 <- brm(formula=nlform, data=carn403, family = gaussian(), prior = nlprior, iter=4000, warmup=2000, chains=4, control = list(adapt_delta = 0.999), seed=15, file=here("fits","fit403"))

# calculating the regression lines for the individual models: 

# 80C

time.seq <- data.frame(time = seq(from = 0, to = 10, by = 0.4))

fitline1 <-
  fitted(model353, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predline1 <-
  predict(model353, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline1 <- ggplot(data = carn353, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline1, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline1, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline1, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "353 K")

#85C

time.seq <- data.frame(time = seq(from = 0, to = 10, by = 0.4))

fitline2 <-
  fitted(model358, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predline2 <-
  predict(model358, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline2 <- ggplot(data = carn358, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline2, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline2, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline2, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "358 K")

#90 C

time.seq <- data.frame(time = seq(from = 0, to = 6, by = 0.2))

fitline3 <-
  fitted(model363, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

predline3 <-
  predict(model363, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline3 <- ggplot(data = carn363, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline3, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline3, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline3, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "363 K")

#95 C

time.seq <- data.frame(time = seq(from = 0, to = 5, by = 0.2))

fitline4 <-
  fitted(model368, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predline4 <-
  predict(model368, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline4 <- ggplot(data = carn368, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline4, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline4, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline4, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "368 K")

# 100 C

time.seq <- data.frame(time = seq(from = 0, to = 4, by = 0.2))

fitline5 <-
  fitted(model373, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predline5 <-
  predict(model373, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline5 <- ggplot(data = carn373, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline5, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline5, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline5, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "373 K")

#110 C

time.seq <- data.frame(time = seq(from = 0, to = 3, by = 0.1))

fitline6 <-
  fitted(model383, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predline6 <-
  predict(model383, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline6 <- ggplot(data = carn383, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline6, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline6, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline6, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "383 K")

#120

time.seq <- data.frame(time = seq(from = 0, to = 2, by = 0.05))

fitline7 <-
  fitted(model393, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predline7 <-
  predict(model393, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline7 <- ggplot(data = carn393, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline7, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline7, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline7, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "393 C")

# 130 C
time.seq <- data.frame(time = seq(from = 0, to = 1, by = 0.05))

fitline8 <-
  fitted(model403, 
         newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

# calculating the prediction lines

predline8 <-
  predict(model403, 
          newdata = time.seq) %>%
  as_tibble() %>%
  bind_cols(time.seq)

regrline8 <- ggplot(data = carn403, 
                    aes(x = time, y = conc)) +
  geom_ribbon(data = predline8, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "lightblue") +
  geom_ribbon(data = fitline8, 
              aes(y = Estimate, ymin = `Q2.5`, ymax = `Q97.5`),
              fill = "blue", alpha=0.3) +
  geom_line(data = fitline8, 
            aes(y = Estimate), size = 1/4) +
  geom_point(color = "navyblue", shape = 19, size = 1.5, alpha = 1) +
  theme_bw()+
  labs(x="time (h)", y="carnitin (mg/L)", subtitle = "403 K")

regrline1 + regrline2 + regrline3 + regrline4 +  regrline5 + regrline6 +regrline7 + regrline8+plot_layout(ncol=3)

# code for multilevel trace plots:
model_all_multi_post <- posterior_samples(model_all_multi) %>% select(-lp__)
my_labeller <- as_labeller(x=vars(), default = label_parsed)

p1 <- bayesplot::mcmc_trace(model_all_multi_post, facet_args = list(labeller=my_labeller))+
  theme(axis.text.x = element_text(angle=70, hjust=1))

myfacets <-bayesplot:: facet_bg(fill = "gray", color = NA) +
  bayesplot::facet_text(color = "black", size = 8)

p1+myfacets

# code for plotting group level and population level

# combines the data with predictions using the random effects with re_formula = NULL
newvary1 <- expand.grid(time=seq(from = 0, to =10, by=0.5),temp=c(353, 358))
newvary2 <- expand.grid(time=seq(from = 0, to =6, by=0.25),temp=c(363, 368, 373))
newvary3 <- expand.grid(time=seq(from = 0, to =3, by=0.1),temp=c(383, 393))
newvary4 <- expand.grid(time=seq(from = 0, to =1, by=0.05),temp=c(403))
newvary <- rbind(newvary1,newvary2,newvary3, newvary4)

carn_ind <- cbind(newvary, predict(model_all_multi, newdata=newvary, re_formula = NA)[,-2]) # fit for the population level
carn_ind2 <- cbind(newvary, predict(model_all_multi, newdata=newvary, re_formula = NULL)[,-2]) #fit for the group level
names(carn_ind) <- c("time", "temp", "conc", "lower", "upper")

carn_ind$temp=as.factor(carn_ind$temp)

(ind_fit <-  ggplot(carnitin, aes(x=time, y=conc))+
    geom_point(colour = "#2c3e50",
               fill = "#2c3e50")+
    facet_wrap(~temp, ncol=4, scales="free_x")+
    geom_line(data = carn_ind, aes(y = conc), size = 1, colour="blue") +
    geom_line(data = carn_ind2, aes(y = Estimate), size = 1, colour="red") + #fit for the group level
    labs(x="time (h)", y="carnitin (mg/L)"))

# code for the strace plots of the single-level nth-order model:

model_all_nth2_post <- posterior_samples(model_all_nth2) %>% select(-lp__)
bayesplot::mcmc_trace(model_all_nth2_post)

# code for the correlation plot of the single-level nth-order model:
cor_alldata2 <- dplyr::select(model_all_nth2_post,b_c0_Intercept:sigma)
# change the names of the columns to be displayed in the panels
cor_alldata2 <- setNames(cor_alldata2, c(expression(paste("c"[0])), expression(paste("n"[t])), expression(paste("k"[ref])), expression(paste("E"[a],"/RT"[ref])), expression(sigma[e])))

# use ggally for a pairs plot
(cor_alldata2_plot <-cor_alldata2  %>% 
    ggpairs(diag=list(continuous="densityDiag"),
            mapping=aes(fill="red"),
            upper = list(continuous = wrap("cor", size = 4, stars=FALSE, title="corr. coef.")), 
            labeller=label_parsed)+ 
    theme(strip.text.x = element_text(size = 10, color = "red"),
          strip.text.y = element_text(size = 10, color = "red")) +
    theme(axis.text.x = element_text(angle=70, hjust=1))
)


# code for producing papaja table:
p1 <- summary(model_all_nth2)
summary_p1 <- rbind(data.frame(p1$fixed) %>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS), data.frame(p1$spec_pars)%>% dplyr::select(-Rhat, -Bulk_ESS, -Tail_ESS))

# back transform Ea:
summary_p1[4,] <- summary_p1[4,]*8.314*373.5747/1000

rownames(summary_p1) <- c("$c_0 \\text { (mg/dm}^3)$", "$n_t \\text { (-)}$", "$k_{ref} \\text { ((dm}^3 \\text{mg}^{-1})^{n_t-1} \\text {h}^{-1})$","$E_a \\text{ (kJ/mol)}$", "$\\sigma_e$")
colnames(summary_p1) <- c("mean","SE", "lower bound", "upper bound")

apa_table(
  summary_p1,
  placement = "H",
  align = c("c", "c", "c", "c"),
  caption = "(ref:all-data-tab)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2,2,2)
  ),
  escape = FALSE
)

# code for residual plot of multilevel model:

(resplot1 <- carnitin %>% add_residual_draws(model_all_multi) %>% ggplot(aes(x=.row, y=.residual)) + stat_pointinterval()+facet_wrap(~temp, scales="free_x", ncol=4) + geom_hline(yintercept =0))

# code for QQ plot of multilevel model:

(qqplot1 <- carnitin %>%
    add_residual_draws(model_all_multi) %>%
    mean_qi() %>%
    ggplot(aes(sample = .residual)) +
    geom_qq() +
    geom_qq_line()+
    facet_wrap(~temp))

# FREQUENTIST REGRESSION

# code for frequentist nonlinear regression 

# This analysis is based upon the workflow presented in Granville Matheson's blog: https://www.granvillematheson.com/post/nonlinear-modelling-using-nls-nlme-and-brms/

# regression equation:
order_func <- function(time, c0, nt, kr) {
  (c0^(1-nt)+(nt-1)*kr*time)^(1/(1-nt))
}

# this is the general function to do nls regression:
carnitin_nls_fit_func <- function(pf_df) {
  nls_multstart(conc~order_func(time, c0, nt, kr),  data=pf_df,
                lower=c(c0=50, nt=1, kr=0.1),
                upper=c(c0=80,nt=3, kr=10),
                start_lower = c(c0=50, nt=0.5, kr=0.1),
                start_upper = c(c0=70, nt=3, kr=10),
                iter=100,
                supp_errors="Y")
}
# then apply this to the carnitin data, first nest time and concentration in each temperature group, called carnitinnest
carnitinnested <- carnitin %>% group_by(temp) %>% nest(carnitinnest=c(time, conc))

# then call the nls function for each group; the ~ signifies a function defined on the fly, . denotes the current dataframe for that function on the fly:

carnitinnested <- carnitinnested %>%  mutate(carnitin_nls_fit=map(carnitinnest, ~carnitin_nls_fit_func(.x)))

# map has done as many regressions as there are groups,and has added them as a column carnitin_nls_fit in carnitinnested. Here is one individual plot:
plot_nls <- function(nls_object, data){
  predframe <- tibble(time=seq(from=min(data$time), to=max(data$time), length.out=100)) %>% 
    mutate(conc=predict(nls_object, newdata=list(time=.$time)))
  ggplot(data, aes(x=time, y=conc))+
    geom_point(size=1)+
    geom_line(data=predframe, aes(x=time, y=conc))
  
}
# plots of separate regressions
plot_nls1 <- plot_nls(carnitinnested$carnitin_nls_fit[[1]], carnitinnested$carnitinnest[[1]])+labs(x="time (h)", y="carnitin (mg/L)", subtitle = "T=353 K")
plot_nls2 <- plot_nls(carnitinnested$carnitin_nls_fit[[2]], carnitinnested$carnitinnest[[2]])+labs(x="time (h)", y="carnitin (mg/L)", subtitle = "T=358 K")
plot_nls3 <- plot_nls(carnitinnested$carnitin_nls_fit[[3]], carnitinnested$carnitinnest[[3]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=363 K")
plot_nls4 <- plot_nls(carnitinnested$carnitin_nls_fit[[4]], carnitinnested$carnitinnest[[4]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=368 K")
plot_nls5 <- plot_nls(carnitinnested$carnitin_nls_fit[[5]], carnitinnested$carnitinnest[[5]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=373 K")
plot_nls6 <- plot_nls(carnitinnested$carnitin_nls_fit[[6]], carnitinnested$carnitinnest[[6]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=383 K")
plot_nls7 <- plot_nls(carnitinnested$carnitin_nls_fit[[7]], carnitinnested$carnitinnest[[7]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=393 K")
plot_nls8 <- plot_nls(carnitinnested$carnitin_nls_fit[[8]], carnitinnested$carnitinnest[[8]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=403 K")

plot_nls1 + plot_nls2 + plot_nls3 + plot_nls4+plot_nls5 + plot_nls6 + plot_nls7 + plot_nls8 +plot_layout(ncol=3)

# Code for plotting frequentist parameter estimates from individual regressions:

carnitin_nls_outcomes <- carnitinnested %>% 
  mutate(outpars = map(carnitin_nls_fit, ~broom::tidy(.x))) %>% 
  select(-carnitinnest, -carnitin_nls_fit) %>% 
  unnest(cols="outpars")

carnitin_nls_nt <-  dplyr::filter(carnitin_nls_outcomes, term=="nt")
carnitin_nls_c0 <-  dplyr::filter(carnitin_nls_outcomes, term=="c0")
carnitin_nls_kr <-  dplyr::filter(carnitin_nls_outcomes, term=="kr")

nls_nt_values <- carnitin_nls_nt[,4]
nls_nt_se <- carnitin_nls_nt[,5]
nt_range <- data.frame(carnitin_nls_nt$temp, nls_nt_values, nls_nt_values-nls_nt_se, nls_nt_values+nls_nt_se)

nls_c0_values <- carnitin_nls_c0[,4]
nls_c0_se <- carnitin_nls_c0[,5]
c0_range <- data.frame(carnitin_nls_c0$temp, nls_c0_values, nls_c0_values-nls_c0_se, nls_c0_values+nls_c0_se)

nls_kr_values <- carnitin_nls_kr[,4]
nls_kr_se <- carnitin_nls_kr[,5]
kr_range <- data.frame(carnitin_nls_kr$temp, log(nls_kr_values), log(nls_kr_values-nls_kr_se), log(nls_kr_values+nls_kr_se))

nt_range_plot <- ggplot() + 
  geom_errorbarh(data=nt_range, aes(y=carnitin_nls_nt.temp, x=estimate, xmin=estimate.1, xmax=estimate.2), height=0.2, size=1, color="blue")  +
  geom_point(data=nt_range, mapping=aes(y=carnitin_nls_nt.temp, x=estimate), size=4, shape=21, fill="white") +
  labs(x=expression(paste(n[t],"-values")), y="Temperature ")

c0_range_plot <- ggplot() + 
  geom_errorbarh(data=c0_range, aes(y=carnitin_nls_c0.temp, x=estimate, xmin=estimate.1, xmax=estimate.2), height=0.2, size=1, color="blue")  +
  geom_point(data=c0_range, mapping=aes(y=carnitin_nls_c0.temp, x=estimate), size=4, shape=21, fill="white") +
  labs(x=expression(paste(c[0],"-values")), y="Temperature ")

kr_range_plot <- ggplot() + 
  geom_errorbarh(data=kr_range, aes(y=carnitin_nls_kr.temp, x=estimate, xmin=estimate.1, xmax=estimate.2), height=0.2, size=1, color="blue")  +
  geom_point(data=kr_range, mapping=aes(y=carnitin_nls_kr.temp, x=estimate), size=4, shape=21, fill="white") +
  labs(x=expression(paste("ln(k"[r],")-values")), y="Temperature ")

c0_range_plot+nt_range_plot+kr_range_plot

# Multilevel modeling the frequentist way using nlme:

order_func2 <- function(time, temp, c0, nt, kref, Ea) {
  (c0^(1-nt)+(nt-1)*kref*exp(Ea*(1-373.5747/temp))*time)^(1/(1-nt))
}

newtimes1 <- expand.grid(time=seq(from = 0, to =10, by=0.5),temp=c(353, 358))
newtimes2 <- expand.grid(time=seq(from = 0, to =6, by=0.25),temp=c(363, 368, 373))
newtimes3 <- expand.grid(time=seq(from = 0, to =3, by=0.1),temp=c(383, 393))
newtimes4 <- expand.grid(time=seq(from = 0, to =1, by=0.05),temp=c(403))
carn_predtimes <- rbind(newtimes1,newtimes2,newtimes3, newtimes4)

carnitin_nlme_fit <- nlme(conc ~ order_func2(time, temp, c0, nt, kref, Ea), 
                          data = carnitin,
                          fixed=c0 + nt + kref + Ea ~ 1, 
                          random = nt ~ 1, 
                          groups = ~ temp, 
                          start = c(c0=66.0, nt=1.1, kref=0.4, Ea=27.0),
                          verbose = F)

carn_nlmepreds <- carn_predtimes %>% mutate(.fitted=predict(carnitin_nlme_fit, newdata=carn_predtimes))

ggplot(carnitin, aes(x=time, y=conc))+
  geom_point()+
  geom_line(data=carn_nlmepreds, aes(y=.fitted))+
  facet_wrap(~temp, ncol=4, scales="free_x")


# Putting the nlme results in a table:

nlme_sum1 <- as.data.frame(coef(summary(carnitin_nlme_fit)))

# to extract random coefficients, we can use the function VarCorr which gives the covariance matrix of the random effects:
nlme_sum2 <- VarCorr(carnitin_nlme_fit)
nlme_sum3 <- summary(carnitin_nlme_fit)
nlme_re1 <- c(nlme_sum1[1,1], nlme_sum1[2,1], nlme_sum1[3,1], nlme_sum1[4,1]*8.3147*373.5747/1000,nlme_sum2[1,2], nlme_sum2[2,2])
nlme_re2 <- c(nlme_sum1[1,2],nlme_sum1[2,2], nlme_sum1[3,2], nlme_sum1[4,2]*8.3147*373.5747/1000, NA, NA)
nlme_re1 <- as.numeric(nlme_re1)
nlme_re <- data.frame(nlme_re1, nlme_re2)

rownames(nlme_re) <- c("$c_0 \\text { (g/dm}^3)$", "$n_t \\text { (-)}$", "$k_{ref} \\text { (dm}^3 \\text{mol}^{-1})^{n_t-1} \\text {h}^{-1}$","$E_a \\text{ (kJ/mol)}$", "$\\sigma_{n_t}$","$\\sigma_e$")
colnames(nlme_re) <- c("mean","SE")

apa_table(
  nlme_re,
  placement = "H",
  align = c("c", "c", "c"),
  caption = "(ref:nlme-summarytable)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)

# Classical two-step analysis in the frequentist way:

# use of purrr::map to linear regression in one go:

# transform conc in log(conc):
ln_fun <- function(x) {
  return(transform(x,log(x)))
} 

carn_ln <- carnitin %>% transform(conc=log(conc))
# group data by temperature
carn_ln_by_T <- carn_ln %>% group_by(temp)
# put the grouped data in a list via nest()
nested_car_ln <- carn_ln_by_T %>% nest()

# define the regression function
lm_fun <- function(data) lm(conc~1+time, data=data)

lm_carn_ln <- nested_car_ln  %>% mutate(model=map(data,lm_fun))

carn_ln_outcomes <- lm_carn_ln %>% mutate(outpars=map(model,~broom::tidy(.x))) %>% unnest(cols="outpars") %>% select(-data, -model,-statistic,-p.value) %>% mutate(estimate=-estimate)

carn_ln_outcomes <- carn_ln_outcomes %>% filter(term=="time") 

# Putting the resultsin a table:
Tinv <- c(1/(8.314*353),1/(8.314*358),1/(8.314*363),1/(8.314*368),1/(8.314*373),1/(8.314*383),1/(8.314*393),1/(8.314*403))
# linear estimates of rate constants are in carn_ln_outcomes
lnkr_lm <- log(carn_ln_outcomes$estimate)
lm_df <- data.frame(Tinv,lnkr_lm)

lm_result <- lm(lnkr_lm~1+Tinv, data=lm_df)
lm_estimates <- broom::tidy(lm_result) %>% select(-statistic,-p.value)

lm_estimates[2,2] <- (lm_estimates[2,2])/-1000
lm_estimates[2,3] <- (lm_estimates[2,3])/1000
lm_estimates <- lm_estimates %>% select(-term)

rownames(lm_estimates) <- c("$\\ln k_0 \\text { (h}^{-1})$","$E_a \\text { (kJ/mol)}$")
colnames(lm_estimates) <- c("estimate","SE")

apa_table(
  lm_estimates,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:linArrh-table)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)

# Plotting the Arrhenius results the frequentist way:

lm_df_plot <- lm_df %>% ggplot(aes(x=Tinv, y=lnkr_lm))+
  geom_point()+
  stat_smooth(method = "lm", col = "red")+
  labs(x="1/T", y=expression(paste("ln(k[r])")))

lm_resid_plot <- ggResidpanel::resid_panel(lm_result, plots=c("resid", "qq"))

lm_df_plot / lm_resid_plot

# nonlinear least squares regression to obtain rate constants without log transformation and plotting the results:

# regression equation: for first-order model
order_func3 <- function(time, c0, kr) {
  (c0*exp(-kr*time))
}

# this is the general function to do nls regression:
carnitin_fo_fit_func <- function(pf_df) {
  nls_multstart(conc~order_func3(time, c0, kr),  data=pf_df,
                lower=c(c0=50, kr=0.1),
                upper=c(c0=80,kr=10),
                start_lower = c(c0=50, kr=0.1),
                start_upper = c(c0=70, kr=10),
                iter=100,
                supp_errors="Y")
}
# first nest time and concentration in each temperature group, call it carnitinnest
carnitin_fo_nested <- carnitin %>% group_by(temp) %>% nest(carnitin_fo_nest=c(time, conc))

# then call the nls function for each group; the ~ signifies a function defined on the fly, . denotes the current dataframe for that function on the fly:

carnitin_fo_nested <- carnitin_fo_nested %>%  mutate(carnitin_fo_fit=map(carnitin_fo_nest, ~carnitin_fo_fit_func(.x)))

# map has done as many regressions as there are groups,and has added them as a column carn_nls_fit in carnitinnested. Here is one individual plot:
plot_fo <- function(fo_object, data){
  predframe <- tibble(time=seq(from=min(data$time), to=max(data$time), length.out=100)) %>% 
    mutate(conc=predict(fo_object, newdata=list(time=.$time)))
  ggplot(data, aes(x=time, y=conc))+
    geom_point(size=1)+
    geom_line(data=predframe, aes(x=time, y=conc))
  
}
# plots of separate regressions
plot_fo1 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[1]], carnitin_fo_nested$carnitin_fo_nest[[1]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=353 K")
plot_fo2 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[2]], carnitin_fo_nested$carnitin_fo_nest[[2]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=358 K")
plot_fo3 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[3]], carnitin_fo_nested$carnitin_fo_nest[[3]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=363 K")
plot_fo4 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[4]], carnitin_fo_nested$carnitin_fo_nest[[4]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=368 K")
plot_fo5 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[5]], carnitin_fo_nested$carnitin_fo_nest[[5]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=373 K")
plot_fo6 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[6]], carnitin_fo_nested$carnitin_fo_nest[[6]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=383 K")
plot_fo7 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[7]], carnitin_fo_nested$carnitin_fo_nest[[7]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=393 K")
plot_fo8 <- plot_fo(carnitin_fo_nested$carnitin_fo_fit[[8]], carnitin_fo_nested$carnitin_fo_nest[[8]])+labs(x="time (h)", y="carnitin (mg/L)",subtitle = "T=403 K")

plot_fo1 + plot_fo2 + plot_fo3 + plot_fo4+plot_fo5 + plot_fo6 + plot_fo7 + plot_fo8 +plot_layout(ncol=3)

# collecting the results and plotting:
carnitin_fo_outcomes <- carnitin_fo_nested %>% 
  mutate(outpars = map(carnitin_fo_fit, ~broom::tidy(.x))) %>% 
  select(-carnitin_fo_nest, -carnitin_fo_fit) %>% 
  unnest(cols="outpars")

carnitin_fo_c0 <-  dplyr::filter(carnitin_fo_outcomes, term=="c0")
carnitin_fo_kr <-  dplyr::filter(carnitin_fo_outcomes, term=="kr")

fo_c0_values <- carnitin_fo_c0[,4]
fo_c0_se <- carnitin_fo_c0[,5]
fo_c0_range <- data.frame(carnitin_fo_c0$temp, fo_c0_values, fo_c0_values-fo_c0_se, fo_c0_values+fo_c0_se)

fo_kr_values <- carnitin_fo_kr[,4]
fo_kr_se <- carnitin_fo_kr[,5]
fo_kr_range <- data.frame(carnitin_fo_kr$temp, log(fo_kr_values), log(fo_kr_values-fo_kr_se), log(fo_kr_values+fo_kr_se))

Tinv <- c(1/(8.314*353),1/(8.314*358),1/(8.314*363),1/(8.314*368),1/(8.314*373),1/(8.314*383),1/(8.314*393),1/(8.314*403))
# nonlinear estimates of rate constants are in fo_kr_values
lnkr_nls <- log(fo_kr_values)
nls_df <- data.frame(Tinv,lnkr_nls)

lm_nls_result <- lm(nls_df$estimate~1+Tinv, data=nls_df)
lm_nls_estimates <- broom::tidy(lm_nls_result) %>%  select(-statistic, -p.value)
lm_nls_estimates[2,2] <- lm_nls_estimates[2,2]/-1000
lm_nls_estimates[2,3] <- lm_nls_estimates[2,3]/1000

nls_df_plot <- nls_df %>% ggplot(aes(x=Tinv, y=estimate))+
  geom_point()+
  stat_smooth(method = "lm", col = "red")+
  labs(x="1/T", y=expression(paste("ln(k[r])")))
nls_resid_plot <- ggResidpanel::resid_panel(lm_nls_result, plots=c("resid", "qq"))

nls_df_plot / nls_resid_plot

# results for papaja table:

Tinv <- c(1/(8.314*353),1/(8.314*358),1/(8.314*363),1/(8.314*368),1/(8.314*373),1/(8.314*383),1/(8.314*393),1/(8.314*403))
# nonlinear estimates of rate constants are in fo_kr_values
lnkr_nls <- log(fo_kr_values)
nls_df <- data.frame(Tinv,lnkr_nls)

lm_nls_result <- lm(nls_df$estimate~1+Tinv, data=nls_df)
lm_nls_estimates <- broom::tidy(lm_nls_result) %>%  select(-statistic, -p.value)
lm_nls_estimates[2,2] <- lm_nls_estimates[2,2]/-1000
lm_nls_estimates[2,3] <- lm_nls_estimates[2,3]/1000

rownames(lm_nls_estimates) <- c("$\\ln k_0$","$E_a \\text{ (kJ/mol)}$")
colnames(lm_nls_estimates) <- c("estimate","SE")

apa_table(
  lm_nls_estimates,
  placement = "H",
  align = c("c", "c"),
  caption = "(ref:nls-linArrh)",
  note = NULL,
  small = TRUE,
  format.args=list(
    digits = c(2,2)
  ),
  escape = FALSE
)



