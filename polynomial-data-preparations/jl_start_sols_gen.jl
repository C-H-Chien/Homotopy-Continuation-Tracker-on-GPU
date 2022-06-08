using HomotopyContinuation

# alea6 has 6 unknowns
@polyvar x1 x2 x3 x4 x5 x6

# define start coefficients (this can be randomly generated)
a1 = 0.7634445720786067 + 0.42887027730523997im
a2 = 0.48885530211853623 + 0.11138063973800438im
a3 = 0.23670283099920075 + 0.970170654114934im
a4 = 0.8570180302948744 + 0.2841877488377911im
a5 = 0.08725168922941018 + 0.48373141757727467im
a6 = 0.7920495030655841 + 0.47288052150374615im
a7 = 0.7383868554797113 + 0.13639129461454735im
a8 = 0.06068556281106052 + 0.05722945887600783im
a9 = 0.5998797349794218 + 0.7565102146522324im
a10 = 0.052962650471145345 + 0.27552639402833656im
a11 = 0.5932690621177386 + 0.25363141425612246im
a12 = 0.24979749791648298 + 0.40070158788540655im
a13 = 0.8967808134473649 + 0.6893566861413845im
a14 = 0.4829398559218605 + 0.16271101220093076im
a15 = 0.016949150420769055 + 0.3917068613819905im
a16 = 0.45377387616324993 + 0.9464622376853775im
a17 = 0.2890666035097773 + 0.6850699215236049im
a18 = 0.013967048107511815 + 0.0037268785987412123im
a19 = 0.821702687633366 + 0.9069173712769043im
a20 = 0.14996850308134202 + 0.34909129115183846im
a21 = 0.6502325834157312 + 0.35558247799378684im
a22 = 0.21046457020819465 + 0.8158953454012239im
a23 = 0.29231860252263697 + 0.12474977691430111im
a24 = 0.4574610057623556 + 0.6637273248038882im
a25 = 0.08454603006957595 + 0.3799252264569122im
a26 = 0.9875522692264518 + 0.9167238406804208im
a27 = 0.8248648454990604 + 0.03507497219930711im
a28 = 0.37248611402298204 + 0.6351619527093473im
a29 = 0.4993021639734443 + 0.32378615791106435im

# start system formulation
f1 = a1*x1^2*x4 + a2*x2*x4*x5 + a3*x2*x4*x6 + a4*x4*x6 + a5*x5*x6
f2 = a6*x1*x2*x6 + a7*x2^2*x5 + a8*x2*x3*x5 + a9*x2*x5^2 + a10*x3^2 + a11*x4*x5*x6
f3 = a12*x1^2*x4 + a13*x1^2 + a14*x1*x4*x6 + a13*x2^2*x5 + a3*x2*x4^2 + a15*x6^3
f4 = a16*x1*x4^2 + a14*x2*x4 + a17*x2*x5^2 + a18*x3*x4 + a19*x3*x5*x6 + a20*x4*x5^2
f5 = a21*x1^2*x3 + a22*x1*x2*x3 + a23*x1*x4*x5 + a24*x1*x4 + a25*x2^2*x5
f6 = a26*x1*x3*x4 + a19*x3^2*x4 + a27*x3^2*x6 + a28*x3 + a29*x4^3 + a25*x5

F = [f1, f2, f3, f4, f5, f6]

# solve the start system
result = solve(F)

# write the resultant start solutions to a file and store under the user-defined directory
hc_start_sols = solutions(result)
io = open("/your/start-sols/dir/julia-start-sols-raw", "w");
using DelimitedFiles
writedlm(io, hc_start_sols);
close(io)