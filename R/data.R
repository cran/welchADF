##' Number of puzzles 42 students were able to solve
##' 
##' An artificial dataset created by Lix et al. which recreates data reported by a real study
##' on perception and concentration, on which 42 students were given several puzzles to be solved. 
##' The students are divided into three balanced groups as they had previously been asked to 
##' imagine solving puzzles in the distant future, near future, or not to imagine anything at all 
##' (control group). 
##' 
##' @format A data frame with 42 rows and 2 variables:
##' \describe{
##'   \item{Group}{group of the student (what the student was asked to imagine)}
##'   \item{y}{number of puzzles the student was able to solve, out of 12}
##' }
##' @references Forster, J., Liberman, N., & Friedman, R.S. (2004). Temporal construal effects on
##' abstract and concrete thinking: consequences for insight and creative cognition.
##' Journal of Personality and Social Psychology, 87, 2, 177-189.
##' 
##' @source \url{http://supp.apa.org/psycarticles/supplemental/met_13_2_110/Example_1_%20OneWay_Independent_Groups_Design.pdf}
##' 
##' @examples 
##' omnibus_LSM <- welchADF.test(perceptionData, response = "y", between.s = "Group")
##' omnibus_trimmed <- update(omnibus_LSM, trimming = TRUE)
##' pairwise_LSM <- update(omnibus_LSM, effect = "Group", contrast = "all.pairwise")
##' pairwise_trimmed <- update(pairwise_LSM, trimming = TRUE, effect.size = TRUE)
##' pairwise_trimmed_boot <- update(pairwise_trimmed, bootstrap = TRUE, seed = 12345, numsim_b = 600)
##' summary(omnibus_LSM)
##' summary(pairwise_trimmed_boot, digits = 6) # digits defaults to max(4, getOption("digits"))
"perceptionData"

##' Students' scores (men and women) on an arithmetic test
##' 
##' An artificial dataset created by Lix et al. from summary data presented by Wicherts et al. (2005)
##' (see the vignette). These authors examined the effects of stereotype threat on women's mathematics ability. 
##' Originally there were four different tests administered to study participants (arithmetic,
##' number series, word problems, and sums tests). The dataset contains only scores for the arithmetic test
##' because these scores exhibited a greater magnitude of variance heterogeneity than scores for
##' the other tests. It is an unbalanced design with cell sizes ranging from 45 to 50 participants, and a total
##' sample size of 283. 
##' 
##' @format A data frame with 283 rows and 3 variables:
##' \describe{
##'    \item{condition}{test condition (control, nullified, stereotype threat)}
##'    \item{sex}{the individual's sex (male, female)}
##'    \item{y}{score on the arithmetic test, out of 40}
##' }
##' @references J.Wicherts, C. Dolan, and D. Hessen. Stereotype threat and group differences in test performance: 
##' a question of measurement invariance. Journal of Personality and Social Psychology, 89(5):696-716, 2005.
##' 
##' @source \url{http://supp.apa.org/psycarticles/supplemental/met_13_2_110/Example_2_Factorial_Independent_Groups_Design.pdf}
##' 
##' @examples
##' omnibus_LSM <- welchADF.test(womenStereotypeData, response = "y", between.s =
##'   c("condition", "sex"), contrast = "omnibus", effect = "condition")
##' omnibus_trimmed <- update(omnibus_LSM, trimming = TRUE, effect = NULL) # unset value of "effect"
##' pairwise_LSM <- update(omnibus_LSM, contrast = "all.pairwise", effect = c("condition", "sex"))
##' pairwise_trimmed <- update(pairwise_LSM, trimming = TRUE)
##' pairwise_trimmed_boot <- update(pairwise_trimmed, bootstrap = TRUE, seed = 12345)
##' summary(omnibus_LSM)
##' summary(omnibus_trimmed)
##' summary(pairwise_trimmed)
##' summary(pairwise_trimmed_boot)
"womenStereotypeData"

##' Children's reaction times (milliseconds) to stimuli of different nature, arranged with
##' four response columns. 
##' 
##' The data (Keselman et al., 2003) represent the reaction times in milliseconds 
##' of children with attention-deficit hyperactivity (ADHD) and normal children 
##' when they are presented four kinds of inputs: 
##' a target alone or an arrow stimuli incongruent, congruent and neutral to the target. 
##' According to the authors, the dataset was artificially generated from the summary 
##' measures given in the original study by Jonkman et al. (1999), in groups of 20 
##' and 10 children to create an unbalanced design. 
##' 
##' @format A data frame with 30 rows and 5 variables:
##' \describe{
##'    \item{Group}{whether the child has ADHD or is healty (normal)}.
##'    \item{TargetAlone}{Reaction time (milliseconds) to a target alone.}
##'    \item{Congruent}{Reaction time (milliseconds) to a congruent stimulus.}
##'    \item{Neutral}{Reaction time (milliseconds) to a neutral stimulus.}
##'    \item{Incongruent}{Reaction time (milliseconds) to an incongruent stimulus.}
##' }
##' @references L. Jonkman, C. Kemner, M. Verbaten, H. van Engeland, J. Kenemans, G. Camfferman, J. Buitelaar, and
##' H. Koelega. Perceptual and response interference in children with attention-deficit hyperactivity
##' disorder, and the effects of methylphenidate. 36(4):419-429, 1999.
##'
##' @source H. J. Keselman, R. R. Wilcox, and L. M. Lix. A generally robust approach to hypothesis testing in
##' independent and correlated groups designs. Psychophyisiology, 40:586-596, 2003. (Data displayed in page 593).
##' 
##' @examples
##' # Assuming a multivariate response
##' omnibus_LSM_multi <- welchADF.test(adhdData, response = c("TargetAlone", "Congruent",
##'   "Neutral", "Incongruent"), between.s = "Group", within.s = "multivariate", contrast = "omnibus")
##' # The same using the S3 method for class formula
##' omnibus_LSM_multi_form <- welchADF.test(cbind(TargetAlone, Congruent, Neutral, Incongruent)
##' ~ Group, data = adhdData)
##'   
##' # Pairwise comparisons of the implicit within-subjects effect present in the multivariate response
##' pairwise_LSM_multi <- update(omnibus_LSM_multi, contrast = "all.pairwise", effect = "multivariate")
##' summary(omnibus_LSM_multi)
##' summary(pairwise_LSM_multi)
"adhdData"

##' Children's reaction times (milliseconds) to stimuli of different nature, arranged with one single
##' response column and taking the multi-variate response as an explicit within-subjects factor.
##' 
##' Exactly the same data explained in "adhdData" but reshaped as follows.
##' 
##' @format A data frame with 120 rows and 4 variables:
##' \describe{
##'    \item{Group}{whether the child has ADHD or is healty (normal)}.
##'    \item{Stimulus}{the stimulus to which the reaction time in this row corresponds.}
##'    \item{Subject}{an integer ID that identifies the subject to which the reaction time corresponds.}
##'    \item{Milliseconds}{reaction time (milliseconds) to of the aforementioned subject to the aforementioned stimulus.}
##' }
##' 
##' @examples
##' # Omnibus test of a mixed between x within subjects model, 
##' # using trimmed means and Winsorized variances
##' omnibus_trimmed <- welchADF.test(adhdData2, response = "Milliseconds",between.s = "Group", 
##'   within.s = "Stimulus", subject = "Subject", contrast = "omnibus", trimming = TRUE)
##'   
##' # The same using S3 method for class formula. The data can be in separate variables in
##' # the scope of the call, not necessarily in a data.frame
##' millis <- adhdData2$Milliseconds
##' gr <- adhdData2$Group
##' st <- adhdData2$Stimulus
##' sbj <- adhdData2$Subject
##' omnibus_trimmed_formula <- welchADF.test(millis ~ gr*st + (st|sbj), 
##'   contrast = "omnibus", trimming = TRUE)
##' summary(omnibus_trimmed_formula)
##' 
##' # Pairwise contrasts of the effects
##' pairwise_LSM <- welchADF.test(adhdData2, response = "Milliseconds", between.s = "Group", 
##'   within.s = "Stimulus", subject = "Subject", contrast = "all.pairwise", effect = "Stimulus")
##' pairwise_trimmed <- update(pairwise_LSM, trimming = TRUE)
##'   
##' # Bootstrapping to obtain an empirical critical value
##' pairwise_trimmed_boot <- update(pairwise_trimmed, bootstrap =TRUE, seed = 123456)
##' summary(pairwise_trimmed_boot)
"adhdData2"

##' Number of visits and time spent in different tunnels of laboratory mice
##' 
##' Wild strain house mice were, at birth, cross fostered onto house mouse (Mus), 
##' deer mouse (Peromyscus) or rat (Rattus) nursing mothers. Ten days after weaning, 
##' each subject was tested in an apparatus that allowed it to enter four different tunnels: 
##' one scented with clean pine shavings, and the other three tunnels with shavings 
##' bearing the scent of Mus, Peromyscus, or Rattus respectively. 
##' Three variables were measured for each tunnel: the number of visits to the 
##' tunnel during a twenty minute test, the time spent by each subject
##' in each of the four tunnels and the latency to first visit of each tunnel.
##' 
##' @format A data frame with 144 rows and 6 variables:
##' \describe{
##'    \item{Subject}{an ID that identifies the mouse to which the measurements of the row correspond.}
##'    \item{nurs}{type of nursing mother this mouse had at birth.}
##'    \item{tunnel}{type of tunnel to which the measurements of this row correspond.}
##'    \item{visits}{number of visits to the aforementioned tunnel.}
##'    \item{time}{total time spent by the mouse at this tunnel.}
##'    \item{latency}{seconds passed before the mouse first visited this tunnel.}
##' }
##' @references K. L.Wuensch. Fostering house mice onto rats and deer mice: Effects on response to species odors.
##' Animal Learning & Behavior, 20(3):253 - 258, 1992
##'
##' @source \url{http://core.ecu.edu/psyc/wuenschk/SPSS/TUNNEL4b.sav}
##' 
##' @examples
##' omnibus_LSM <- welchADF.test(miceData, response = c("visits", "time", "latency"),
##'   between.s = "nurs", within.s = "tunnel", subject = "Subject", contrast = "omnibus")
##' # Formula interface, using cbind() to specify a multivariate response.
##' omnibus_LSM_formula <- welchADF.test(
##'   cbind(visits, time, latency) ~ nurs*tunnel + (tunnel|Subject), data = miceData) 
##' pairwise_LSM_nurs <- welchADF.test(miceData, response = c("visits", "time",
##'   "latency"), between.s = "nurs", within.s = "tunnel", subject = "Subject",
##'   effect = "nurs", contrast = "all.pairwise")
##' pairwise_LSM_tunnel <- update(pairwise_LSM_nurs, effect = "tunnel")
##' \dontrun{
##'   pairwise_nurs_trimmed_boot <- update(pairwise_LSM_nurs, trimming = TRUE, bootstrap = TRUE)
##'   pairwise_tunnel_trimmed_boot <- update(pairwise_nurs_trimmed_boot, effect = "tunnel")
##'   summary(pairwise_nurs_trimmed_boot)
##' }
##' summary(omnibus_LSM)
##' summary(pairwise_LSM_nurs)
##' summary(pairwise_LSM_tunnel)

"miceData"