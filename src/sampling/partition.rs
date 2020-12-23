use super::space::Space;
pub trait EqualSizedGrid {
    type Sp: Space;

    fn new(cell_size: f32) -> Self;
    fn insert(&mut self, p: &<Self::Sp as Space>::Sample) -> VertexIdx;
    fn neighbors(&self, p: &<Self::Sp as Space>::Sample, samples: &[<Self::Sp as Space>::Sample]) -> Vec<&Vec<VertexIdx>>;
}


use crate::{
    sampling::space::TwoDim,
    coord::Point2
};
pub struct Equal2DSizedGrid<F> {
    num_cell_width: usize,
    num_cell_height: usize,
    num_points_inserted: usize,
    grid: Vec<Option<Vec<VertexIdx>>>,
    cell_size: f32,
    f: std::marker::PhantomData<F>,
}
type VertexIdx = usize;

impl<F> EqualSizedGrid for Equal2DSizedGrid<F>
where F: Fn(&Point2) -> bool {
    type Sp = TwoDim<F>;

    fn new(cell_size: f32) -> Self {
        let num_cell_width = (1.0 / cell_size) as usize + 1;
        let num_cell_height = (1.0 / cell_size) as usize + 1;

        let mut grid: Vec<Option<Vec<VertexIdx>>> = vec![None; num_cell_width * num_cell_height];
        let num_points_inserted = 0;
        Self {
            num_cell_width,
            num_cell_height,
            num_points_inserted,
            grid,
            cell_size,
            f: std::marker::PhantomData
        }
    }

    fn insert(&mut self, p: &Point2) -> VertexIdx {
        let i = (p.x / self.cell_size) as usize;
        let j = (p.y / self.cell_size) as usize;
    
        let idx = j * self.num_cell_width + i;
        let idx_new_point = self.num_points_inserted;
        if let Some(idx_points) = &mut self.grid[idx] {
            idx_points.push(idx_new_point);
        } else {
            self.grid[idx] = Some(vec![idx_new_point]);
        }

        self.num_points_inserted += 1;

        idx_new_point
    }

    fn neighbors(&self, p: &Point2, samples: &[Point2]) -> Vec<&Vec<VertexIdx>> {
        // Get the grid idx where the point is
        let idx_col = (p.x / self.cell_size) as i32;
        let idx_row = (p.y / self.cell_size) as i32;

        let idx_col_max = (idx_col + 2).min(self.num_cell_width as i32 - 1) as usize;
        let idx_col_min = (idx_col - 2).max(0) as usize;
        assert!(idx_col_max > idx_col_min);
        let idx_row_max = (idx_row + 2).min(self.num_cell_height as i32 - 1) as usize;
        let idx_row_min = (idx_row - 2).max(0) as usize;
        assert!(idx_row_max > idx_row_min);

        let it = NeighborsIter {
            g: &self,
            cur_row: idx_row_min,
            cur_col: idx_col_min,

            row_rng: idx_row_min..(idx_row_max+1),
            col_rng: idx_col_min..(idx_col_max+1),

            finished: false,
        };

        it.collect()
    }
}
use std::ops::Range;
struct NeighborsIter<'a, F> {
    g: &'a Equal2DSizedGrid<F>,
    cur_row: usize,
    cur_col: usize,

    row_rng: Range<usize>,
    col_rng: Range<usize>,
    finished: bool,
}

impl<'a, F> Iterator for NeighborsIter<'a, F> {
    type Item = &'a Vec<VertexIdx>;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            None
        } else {
            let mut found = None;
            while !self.finished && found.is_none() {
                let idx_grid = self.cur_row * self.g.num_cell_width + self.cur_col;
                found = self.g.grid[idx_grid].as_ref();

                if self.cur_col == self.col_rng.end - 1 {
                    self.cur_row += 1;
                    self.cur_col = self.col_rng.start;
                } else {
                    self.cur_col += 1;
                }

                self.finished = (self.cur_row == self.row_rng.end - 1) && (self.cur_col == self.col_rng.end - 1);
            }

            found
        }
    }
}