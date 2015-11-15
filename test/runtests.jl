using GainPatterns
using Base.Test

#################################################################
# testing angular_error function
#################################################################
@test angular_error(5, 355) == 10
@test angular_error(355, 5) == 10

@test angular_error(5, -5) == 10
@test angular_error(-5, 5) == 10

@test angular_error(340, 17) == 37
@test angular_error(17, 340) == 37

@test angular_error(-17, 340) == 3
@test angular_error(340, -17) == 3

#################################################################
# testing angular_error_rel function
#################################################################
@test angular_error_rel(5, 355) == -10
@test angular_error_rel(355, 5) == 10

@test angular_error_rel(5, -5) == -10
@test angular_error_rel(-5, 5) == 10

@test angular_error_rel(340, 17) == 37
@test angular_error_rel(17, 340) == -37

@test angular_error_rel(-17, 340) == -3
@test angular_error_rel(340, -17) == 3

v1 = [5, 355, 5, -5, 340, 17, -17, 340]
v2 = [355, 5, -5, 5, 17, 340, 340, -17]
vec_ans = [-10, 10, -10, 10, 37, -37, -3, 3]
@test angular_error_rel(v1,v2) == vec_ans
