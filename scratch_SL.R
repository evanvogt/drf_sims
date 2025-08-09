###################
# title: Custom learners for the superlearner
# date started: 22/05/2025
# date finished:
# author: Ellie Van Vogt
###################
# based on code from Matt Pryce


nlambda_seq = c(50,100,250)
alpha_seq <- c(0.5,1)
usemin_seq <- c(FALSE,TRUE)
para_learners = create.Learner("SL.glmnet",
                               tune = list(nlambda = nlambda_seq,
                                           alpha = alpha_seq,
                                           useMin = usemin_seq))
para_learners
#Random forest - One covariate
mtry_seq1 <-  1
min_node_seq <- c(10,20,50)
rf_learners1 = create.Learner("SL.ranger",
                              tune = list(mtry = mtry_seq1,
                                          min.node.size = min_node_seq))
rf_learners1
mtry_seq6 <-  floor(sqrt(6) * c(0.5, 1))
min_node_seq <- c(10,20,50)
rf_learners6 = create.Learner("SL.ranger", tune = list(mtry = mtry_seq6, min.node.size = min_node_seq))
rf_learners6

# scratch for checking the resources for superlearner runs

predict.SuperLearner()
DR_sub$method$computePred()
worker_ppid <- 3542457
worked_ppid <- 3541929
worker_ppid <- 3541952
worker_ppid <- 	2317178
syrup_565 %>%
  filter(ppid == worker_ppid | pid == worker_ppid) %>%
  ggplot() +
  aes(x = id, y = pct_cpu, group = pid) +
  geom_line() +
  scale_x_continuous(breaks = 1:max(syrup_565$id))



SSRMST::ssrmst()