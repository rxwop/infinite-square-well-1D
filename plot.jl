using Plots, QuadGK, ColorSchemes

dx = 0.1
dt = 1/30
elapse = 10
a = 10
m = 0.2

h = 1.054571817
xran = 0:dx:a
yran = 0:dx:elapse
lbls = ["real" "imaginary"]
clrs = ColorSchemes.GnBu_9



# Vector e.g. [1, 4, 5] or collected Range e.g. collect(1:5)
input = [9,15,30]
states = sort(input)

function constant(n)
    1/sqrt(size(states, 1))
    #2/(n*pi)*(1-cos(n*pi/2))
end

Cₙ = [constant(i) for i in states]

Ψₙ(x) = [sin(pi*x*i/a)*sqrt(2/a) for i in states]
dΨₙ(x) = [cos(pi*x*i/a)*sqrt(2/a)*pi*i/a for i in states]

φₙ(t) = [exp(-im*(i*pi)^2*h*t/(2*m*a^2)) for i in states]

radii = [[0,abs(i)] for i in Cₙ]

Eₙ = [(i*pi*h)^2/2*m*a^2 for i in states]

anim = @animate for ti in 0:dt:elapse

        Ψ(x) = sum(Ψₙ(x).*φₙ(ti).*Cₙ)
        dΨ(x) = sum(dΨₙ(x).*φₙ(ti).*Cₙ)

        ref(x) = real(Ψ(x))
        imf(x) = imag(Ψ(x))
        ρ(x) = abs2(Ψ(x))
        J(x) = real(im*h/(2*m)*(Ψ(x)*conj(dΨ(x))-conj(Ψ(x))*dΨ(x)))
        xΨ(x) = x*ρ(x)
        x²Ψ(x) = x^2*ρ(x)
        expcX = first(quadgk(xΨ, 0, a))
        expcX² = first(quadgk(x²Ψ, 0, a))
        σₓ = sqrt(expcX²-expcX^2)

        theta = [[0,i] for i in angle.(φₙ(ti))]

        p1 = plot(xran, [ref, imf], ylims = [-1,1], color = [clrs[5] clrs[3]], labels = lbls, title = "Ψ(x, t)")
        p2 = plot(xran, ρ, ylims = [0,1], color = clrs[6], labels = false, title = "∣Ψ(x, t)∣²", fill = true)
        p3 = plot(theta, radii, proj =:polar, lims = [0,1], showaxis = false, labels = false, palette = [clrs[4], clrs[5], clrs[6], clrs[7], clrs[8]], title = "φ(t)")
        scatter!(p2, [expcX], [0.5], color = clrs[5], xerror = σₓ, markerstrokecolor=clrs[8], labels = false)
        
        p5 = plot(xran, J, ylims = [-1.5, 1.5], labels = false, title = "J(x,t)", color = clrs[6])
        
        l1 = @layout[a b; c{0.4h} d]

        plot(p1, p2, p3, p5, fontfamily = :serif, layout = l1)
        
end

Eplot = vcat(((states[begin]-1)*pi*h)^2/2*m*a^2, Eₙ, ((states[end]+1)*pi*h)^2/2*m*a^2)
Cplot = vcat(0, abs2.(Cₙ), 0)
Psi(x, t) = abs2(sum(Ψₙ(x).*φₙ(t).*Cₙ))
l2 = @layout[a{0.3w} [b; c]]

p4 = plot(Cplot, Eplot, xlims = [0, 1], yticks = :log10, markershapes = :circle, ylabel = "E", xlabel = "∣C∣²", labels = false, line = clrs[4], markercolors = clrs[4])
hline!(p4, [sum(abs2.(Cₙ).*Eₙ)], color = clrs[9], label = false)
p8 = surface(xran, yran, Psi, camera = [15, 75], seriescolor = cgrad(clrs), ylabel = "t", zlabel = "∣Ψ(x, t)∣²")
p6 = contourf(xran, yran, Psi, seriescolor = cgrad(clrs), ylabel = "t")

stats = plot(p4, p8, p6, layout = l2, fontfamily = :serif)

savefig(stats, "isw1_$input _stats.png")

gif(anim, "isw1_$input.gif", fps = 1/dt)