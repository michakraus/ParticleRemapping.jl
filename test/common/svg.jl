using JSON
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
    @test occursin("b &lt;c&gt;", s) && !occursin("b <c>", s)
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
    # clip ids differ between figures, so that several fit on one page
    ids(svg) = [m.captures[1] for m in eachmatch(r"clipPath id=\"(pr-clip-\d+)\"", svg)]
    @test isempty(intersect(ids(s), ids(s2)))
end

@testset "player" begin
    js = svg_player_json(stages; window = (0, 1, -1, 1), id = "t",
        titles = ["a \"q\" \\ b", "b <c>"],
        targets = [nothing, ([0.2, 0.3, 0.4], [0.5, 0.5, 0.5])],
        overlays = k -> [
            svg_line(0, 0, 1, 0; dash = true), svg_arrow(0, 0, 1, 1),
            svg_rect(0, 0, 0.5, 0.5), svg_text(0.5, 0.5, "t$k")])
    d = JSON.parse(js)
    @test d["version"] == 1
    @test d["id"] == "t"
    @test d["window"] == [0.0, 1.0, -1.0, 1.0]
    @test d["markers"] == 3
    steps = d["steps"]
    @test length(steps) == 2

    s1, s2 = steps[1], steps[2]
    @test s1["caption"] == "a \"q\" \\ b"
    @test s2["caption"] == "b <c>"
    @test s1["x"] == [0.1, 0.5, 0.9]
    @test s1["v"] == [-0.5, 0.0, 0.5]
    @test s2["x"] == [0.2, 0.6, 0.8]
    @test s2["v"] == [0.5, 0.1, -0.5]
    for st in steps
        @test length(st["r"]) == 3
        @test length(st["c"]) == 3
        @test all(c -> occursin(r"^#[0-9a-f]{6}\z", c), st["c"])
    end
    # shared colour range and size reference across all states
    @test s1["c"][1] == s2["c"][3]
    @test s1["r"][1] == s2["r"][3]
    @test s1["c"][1] != s2["c"][1]
    @test s1["r"][1] != s2["r"][1]
    @test s1["r"][3] == 7.0
    @test s2["r"][1] == 7.0

    # targets
    @test s1["targets"] === nothing
    @test s2["targets"]["x"] == [0.2, 0.3, 0.4]
    @test s2["targets"]["v"] == [0.5, 0.5, 0.5]

    # all four overlay kinds and per-step evaluation
    for (k, st) in enumerate(steps)
        ovs = st["overlays"]
        @test length(ovs) == 4
        @test ovs[1]["kind"] == "line"
        @test ovs[1]["x1"] == 0.0 && ovs[1]["v1"] == 0.0 &&
              ovs[1]["x2"] == 1.0 && ovs[1]["v2"] == 0.0
        @test ovs[1]["color"] == "#404040" && ovs[1]["width"] == 1.0 &&
              ovs[1]["dash"] === true
        @test ovs[2]["kind"] == "arrow"
        @test ovs[2]["x1"] == 0.0 && ovs[2]["v1"] == 0.0 &&
              ovs[2]["x2"] == 1.0 && ovs[2]["v2"] == 1.0
        @test ovs[2]["color"] == "#808080" && ovs[2]["width"] == 1.0
        @test ovs[3]["kind"] == "rect"
        @test ovs[3]["x"] == 0.0 && ovs[3]["v"] == 0.0 &&
              ovs[3]["dx"] == 0.5 && ovs[3]["dv"] == 0.5
        @test ovs[3]["color"] == "#56B4E9" && ovs[3]["opacity"] == 0.15
        @test ovs[4]["kind"] == "text"
        @test ovs[4]["x"] == 0.5 && ovs[4]["v"] == 0.5 &&
              ovs[4]["text"] == "t$k" && ovs[4]["color"] == "#404040" &&
              ovs[4]["size"] == 11.0
    end

    # category colours, vector overlays, tuple targets, function colours and control-char escaping
    d_cat = JSON.parse(svg_player_json(stages; window = (0, 1, -1, 1), id = "cat",
        colors = svg_category_colors([1, 2, 8]),
        targets = ([0.1, 0.5, 0.9], [0.0, 0.0, 0.0]),
        overlays = [[svg_line(0, 0, 1, 0; dash = false)], []],
        titles = ["line1\n\t\r\b\f\0", "step 2"]))
    @test all(st -> st["c"] == ["#0072b2", "#e69f00", "#0072b2"], d_cat["steps"])
    @test d_cat["steps"][1]["targets"]["x"] == [0.1, 0.5, 0.9]
    @test d_cat["steps"][2]["targets"]["x"] == [0.1, 0.5, 0.9]
    @test d_cat["steps"][1]["overlays"][1]["dash"] === false
    @test isempty(d_cat["steps"][2]["overlays"])
    @test d_cat["steps"][1]["caption"] == "line1\n\t\r\b\f\0"

    # function colours and rgb(r,g,b) normalisation
    d_rgb = JSON.parse(svg_player_json(stages; window = (0, 1, -1, 1), id = "rgb",
        colors = k -> k == 1 ? ["rgb(0,15,255)", "#ABCDEF", "rgb(1, 2, 3)"] : nothing))
    @test d_rgb["steps"][1]["c"] == ["#000fff", "#abcdef", "#010203"]
    @test d_rgb["steps"][2]["c"] != d_rgb["steps"][1]["c"]

    # one state produces one step
    d_one = JSON.parse(svg_player_json([stages[1]]; window = (0, 1, -1, 1), id = "one"))
    @test length(d_one["steps"]) == 1
    @test d_one["steps"][1]["caption"] == "a"

    # asymmetric coordinates on every overlay kind, id and overlay-text escaping, and 4-decimal rounding
    stages_asym = [
        RemapStage(mk([1 / 3, 0.5, 2 / 3], [-0.5, 0.0, 0.5], [0.1, 0.2, 0.3]), "a"),
        RemapStage(mk([0.2, 0.6, 0.8], [0.5, 0.1, -0.5], [0.2, 0.4, 0.6]), "b")]
    d_asym = JSON.parse(svg_player_json(stages_asym; window = (0, 1, -1, 1),
        id = "id \"q\" \\ \x01\n\t",
        overlays = k -> [
            svg_line(0.1, -0.2, 0.8, 0.3; color = "#345678", width = 2.5, dash = true),
            svg_arrow(0.15, -0.25, 0.75, 0.45; color = "#456789", width = 1.5),
            svg_rect(0.12, -0.35, 0.33, 0.67; color = "#123456", opacity = 0.25),
            svg_text(0.22, 0.73, "t$k \"q\" \\ \x02\n\r\b\f"; color = "#234567", size = 13)]))
    @test d_asym["id"] == "id \"q\" \\ \x01\n\t"
    # 4-decimal rounding
    @test d_asym["steps"][1]["x"] == [0.3333, 0.5, 0.6667]
    # shared colour range and size reference with unequal per-state weight extrema/maxima
    a1, a2 = d_asym["steps"][1], d_asym["steps"][2]
    @test a1["c"][2] == a2["c"][1]
    @test a1["r"][2] == a2["r"][1]
    @test a1["r"][3] < 7.0
    @test a2["r"][3] == 7.0
    @test a1["c"][3] != a2["c"][3]
    # every overlay coordinate field and escaped overlay text
    for (k, st) in enumerate(d_asym["steps"])
        ovs = st["overlays"]
        @test ovs[1]["x1"] == 0.1 && ovs[1]["v1"] == -0.2 &&
              ovs[1]["x2"] == 0.8 && ovs[1]["v2"] == 0.3 &&
              ovs[1]["color"] == "#345678" && ovs[1]["width"] == 2.5 &&
              ovs[1]["dash"] === true
        @test ovs[2]["x1"] == 0.15 && ovs[2]["v1"] == -0.25 &&
              ovs[2]["x2"] == 0.75 && ovs[2]["v2"] == 0.45 &&
              ovs[2]["color"] == "#456789" && ovs[2]["width"] == 1.5
        @test ovs[3]["x"] == 0.12 && ovs[3]["v"] == -0.35 &&
              ovs[3]["dx"] == 0.33 && ovs[3]["dv"] == 0.67 &&
              ovs[3]["color"] == "#123456" && ovs[3]["opacity"] == 0.25
        @test ovs[4]["x"] == 0.22 && ovs[4]["v"] == 0.73 &&
              ovs[4]["text"] == "t$k \"q\" \\ \x02\n\r\b\f" &&
              ovs[4]["color"] == "#234567" && ovs[4]["size"] == 13.0
    end

    # refusals (every guard with an anchored message pattern)
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 2$" svg_player_json(
        [
            stages[1], RemapStage(mk([NaN, 0.5, 0.8], [0.5, 0.1, -0.5], [0.3, 0.2, 0.1]), "b")];
        window = (0, 1, -1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 2$" svg_player_json(
        [
            stages[1], RemapStage(mk([0.2, 0.5, 0.8], [0.5, Inf, -0.5], [0.3, 0.2, 0.1]), "b")];
        window = (0, 1, -1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 1$" svg_player_json(
        [
            RemapStage(mk([0.1, 0.5, 0.9], [-0.5, 0.0, 0.5], [Inf, 0.2, 0.3]), "a"), stages[2]];
        window = (0, 1, -1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 2$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        targets = [nothing, ([0.2, 0.3, 0.4], [0.5, NaN, 0.5])])
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 2$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        targets = [nothing, ([Inf, 0.3, 0.4], [0.5, 0.5, 0.5])])
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 1$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        overlays = [[svg_line(NaN, 0, 1, 0)], []])
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 1$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        overlays = [[svg_rect(0, 0, Inf, 0.5)], []])
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 1$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        overlays = [[svg_text(0, 0, "x"; size = NaN)], []])
    @test_throws r"^ArgumentError: svg_player_json: a non-finite number in state 1$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        overlays = [[svg_arrow(0, 0, 1, 1; width = Inf)], []])
    @test_throws r"^ArgumentError: svg_player_json: the window must hold finite bounds with xmin < xmax and vmin < vmax$" svg_player_json(
        stages; window = (0, 0, -1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: the window must hold finite bounds with xmin < xmax and vmin < vmax$" svg_player_json(
        stages; window = (0, 1, 1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: the window must hold finite bounds with xmin < xmax and vmin < vmax$" svg_player_json(
        stages; window = (0, 1, 1, -1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: the window must hold finite bounds with xmin < xmax and vmin < vmax$" svg_player_json(
        stages; window = (0, Inf, -1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: no states$" svg_player_json(
        RemapStage[]; window = (0, 1, -1, 1), id = "t")
    @test_throws r"^DimensionMismatch: svg_player_json: every state must hold the same markers$" svg_player_json(
        [stages[1], RemapStage(mk([0.1], [0.1], [0.1]), "c")];
        window = (0, 1, -1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: the states hold no markers$" svg_player_json(
        [RemapStage(mk(Float64[], Float64[], Float64[]), "empty")];
        window = (0, 1, -1, 1), id = "t")
    @test_throws r"^ArgumentError: svg_player_json: the caption of state 2 is empty$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", titles = ["a", ""])
    @test_throws r"^ArgumentError: svg_player_json: the caption of state 2 is empty$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", titles = ["a", "   "])
    @test_throws r"^DimensionMismatch: svg_player_json: the targets of state 2 must hold one point per marker$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        targets = [nothing, ([0.2, 0.3], [0.5, 0.5])])
    @test_throws r"^ArgumentError: svg_player_json: the colour \"white\" of state 1 is not #rrggbb or rgb\(r,g,b\)$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        colors = ["white", "#000000", "#000000"])
    @test_throws r"^ArgumentError: svg_player_json: the colour \"rgb\(256,0,0\)\" of state 1 is not #rrggbb or rgb\(r,g,b\)$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        colors = ["rgb(256,0,0)", "#000000", "#000000"])
    @test_throws r"^ArgumentError: svg_player_json: the colour \"#123456\\n\" of state 1 is not #rrggbb or rgb\(r,g,b\)$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        colors = ["#123456\n", "#000000", "#000000"])
    @test_throws r"^ArgumentError: svg_player_json: the colour \"rgb\(1,2,3\)\\n\" of state 1 is not #rrggbb or rgb\(r,g,b\)$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        colors = ["rgb(1,2,3)\n", "#000000", "#000000"])
    @test_throws r"^ArgumentError: svg_player_json: an overlay of state 1 has the unknown kind :circle$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t",
        overlays = [[(; kind = :circle)], []])
    @test_throws r"^ArgumentError: svg_player_json: rmax must be finite and positive" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", rmax = 0)
    @test_throws r"^ArgumentError: svg_player_json: rmax must be finite and positive" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", rmax = -1)
    @test_throws r"^ArgumentError: svg_player_json: rmax must be finite and positive" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", rmax = Inf)
    @test_throws r"^ArgumentError: svg_player_json: rmax must be finite and positive" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", rmax = NaN)
    @test_throws r"^DimensionMismatch: svg_player_json: titles must hold one caption per state$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", titles = ["only one"])
    @test_throws r"^DimensionMismatch: svg_player_json: targets must hold one entry per state$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", targets = [nothing])
    @test_throws r"^DimensionMismatch: svg_player_json: overlays must hold one entry per state$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", overlays = [[]])
    @test_throws r"^DimensionMismatch: svg_player_json: the colours of state 1 must hold one colour per marker$" svg_player_json(
        stages; window = (0, 1, -1, 1), id = "t", colors = ["#000000"])
end
