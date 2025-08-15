struct MarchingSquare {}

fn extract_borders<F>(num_sampling_vertices: usize, inside_area: F)
where
    F: Fn(Point2) -> bool,
{
    let square_size = 1.0 / (num_sampling_vertices as f32 - 1.0);

    let mut vertices = vec![];
    let mut segments: Vec<(usize, usize)> = vec![];
    let mut g = HashMap::new();
    for i in 0..(num_sampling_vertices - 1) {
        let y = (i as f32) * square_size;
        for j in 0..(num_sampling_vertices - 1) {
            let x = (j as f32) * square_size;

            let tl_in = inside_area(Point2::new(x, y));
            let tr_in = inside_area(Point2::new(x + square_size, y));
            let bl_in = inside_area(Point2::new(x, y + square_size));
            let br_in = inside_area(Point2::new(x + square_size, y + square_size));

            let id_segment = segments.len();
            let id_vertex = vertices.len();

            match (tl_in, tr_in, bl_in, br_in) {
                (true, true, true, true) => (),
                (false, false, false, true) => {
                    vertices.push(Point2::new(x + square_size, y + 0.5 * square_size));
                    vertices.push(Point2::new(x + 0.5 * square_size, y + square_size));

                    g.get_mut(&(i, j)).push(id_segment)
                    segments.push((id_vertex, id_vertex + 1));
                },
                (false, false, true, false) => {
                    vertices.push(Point2::new(x, y + 0.5 * square_size));
                    vertices.push(Point2::new(x + 0.5 * square_size, y + square_size));

                    g.get_mut(&(i, j)).push(id_segment)
                    segments.push((id_vertex, id_vertex + 1));
                },
                (false, true, false, false) => {
                    vertices.push(Point2::new(x + 0.5 * square_size, y));
                    vertices.push(Point2::new(x + square_size, y + 0.5*square_size));

                    g.get_mut(&(i, j)).push(id_segment)
                    segments.push((id_vertex, id_vertex + 1));
                },
                (true, false, false, false) => {
                    vertices.push(Point2::new(x + 0.5 * square_size, y));
                    vertices.push(Point2::new(x, y + 0.5*square_size));

                    g.get_mut(&(i, j)).push(id_segment)
                    segments.push((id_vertex, id_vertex + 1));
                },
                (false, false, true, true) => {
                    vertices.push(Point2::new(x, y + 0.5*square_size));
                    vertices.push(Point2::new(x + square_size, y + 0.5*square_size));

                    g.get_mut(&(i, j)).push(id_segment)
                    segments.push((id_vertex, id_vertex + 1));
                },
                _ => ()
            }
        }
    }
    
}
