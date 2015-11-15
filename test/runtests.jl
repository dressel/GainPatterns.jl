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
