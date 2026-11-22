using ParticleRemapping
using Test

mk(x, v, w) = Markers(x, v, w; m = 1.0)
stages = [RemapStage(mk([0.1, 0.5, 0.9], [-0.5, 0.0, 0.5], [0.1, 0.2, 0.3]), "a"),
    RemapStage(mk([0.2, 0.6, 0.8], [0.5, 0.1, -0.5], [0.3, 0.2, 0.1]), "b <c>")]
count_of(s, pattern) = length(collect(eachmatch(pattern, s)))

@testset "panels" begin
    f = svg_stage_panels(stages; window = (0, 1, -1, 1), columns = 2,
        overlays = k -> [
            svg_line(0, 0, 1, 0), svg_arrow(0, 0, 1, 1), svg_rect(0, 0, 0.5, 0.5),
            svg_text(0.5, 0.5, "t$k")],
        targets = [nothing, ([0.2], [0.5])])
    s = sprint(show, MIME("image/svg+xml"), f)
    @test startswith(s, "<svg") && endswith(s, "</svg>")
    @test occursin("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>", s)
    @test count_of(s, r"<circle") == 6 + 1          # three markers per panel and one target
    @test count_of(s, r"<clipPath id=\"pr-clip-\d+\"") == 2
    @test count_of(s, r"marker-end=") == 2          # one arrow per panel
    @test occursin("t1", s) && occursin("t2", s)
    # a marker at x = 0.1 of a box of width 220 px from x0 = 24 px
    @test occursin("cx=\"46.0\"", s)
    # category colours replace the weight colours
    c = svg_category_colors([1, 2, 8])
    @test c == ["#0072B2", "#E69F00", "#0072B2"]
    s2 = svg_stage_panels(stages; window = (0, 1, -1, 1), colors = c).svg
    @test occursin("fill=\"#E69F00\"", s2) && !occursin("rgb(", s2)
end

@testset "animation" begin
    a = svg_stage_animation(stages; window = (0, 1, -1, 1), hold = 1.0, move = 1.5).svg
    @test count_of(a, r"<circle") == 3
    @test occursin("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>", a)
    @test occursin("b &lt;c&gt;", a) && !occursin("b <c>", a)
    # four animated attributes per marker, two states and a return to the first
    @test count_of(a, r"attributeName=\"cx\"") == 3
    @test occursin("dur=\"5.0s\"", a)
    @test occursin("keyTimes=\"0.0;0.2;0.5;0.7;1.0\"", a)
    # the second title is hidden without animation
    @test occursin("<g opacity=\"0\">", a)
    @test_throws DimensionMismatch svg_stage_animation(
        [stages[1], RemapStage(mk([0.1], [0.1], [0.1]), "c")]; window = (0, 1, -1, 1))
    # clip ids differ between figures, so that several fit on one page
    ids(s) = [m.captures[1] for m in eachmatch(r"clipPath id=\"(pr-clip-\d+)\"", s)]
    @test isempty(intersect(ids(a), ids(svg_stage_animation(stages;
        window = (0, 1, -1, 1)).svg)))
end
