using Gaston, SpecialFunctions
x = y = 0:0.075:10
surf(x, y, (x,y) -> besselj0(y)*x^2, with = "pm3d",
     Axes(view = (45, 45),
          pm3d = "lighting primary 0.5 specular 0.4",
          key = :off)
     )
x = 1:0.1:10
plot(x, x.^0.5,
     w = "l",
     legend = "'Pow 0.5'",
     dt = 2,
     lw = 2,
     lc = :red,
     Axes(grid = :on,
          key = "left",
          axis = "semilogy"))
plot!(x, x,
      w = :l,
      leg = :Pow_1,
      dt = 1,
      lw = 3,
      lc = :blue)
plot!(x, x.^2,
      curveconf = "w l tit 'Pow 2' dt 3 lw 2 lc 'purple'")
