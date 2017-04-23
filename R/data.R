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
##' @source \url{http://homepage.usask.ca/~lml321/Example1.pdf}
"perceptionData"

##' Students' scores (men and women) on an arithmetic test
##' 
##' An artificial dataset created by Lix et al. from summary data presented by Wicherts et al. (2005)
##' (see the vignette). These authors examined the effects of stereotype threat on women's mathematics ability. 
##' Originally there were four different tests administered to study participants (arithmetic,
##' number series, word problems, and sums tests); the dataset contains only scores for the arithmetic test
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
##' @source \url{http://homepage.usask.ca/~lml321/Example2.pdf}
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
"miceData"