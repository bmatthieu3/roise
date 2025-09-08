use crate::Point2;
use std::hash::Hash;

struct Aabb {
    xy_min: Point2<f32>,
    xy_max: Point2<f32>,
}

trait Spatial {
    type Id: Copy + Eq + Hash;

    /// Bounding box (for fast insertion/search)
    fn aabb(&self) -> Aabb;

    /// Optional: exact geometry check if needed
    fn contains(&self, point: Point2<f32>) -> bool;

    fn id(&self) -> Self::Id;
}

use std::collections::HashMap;
struct SpatialGrid<T>
where
    T: Spatial,
{
    cell_size: f32,
    cells: HashMap<(i32, i32), Vec<T::Id>>,

    objects: HashMap<T::Id, T>,
}

impl<T> SpatialGrid<T>
where
    T: Spatial,
{
    fn insert(&mut self, obj: T) {
        let Aabb { xy_min, xy_max } = obj.aabb();
        let id = obj.id();

        let i_min = (xy_min.x / self.cell_size) as i32;
        let j_min = (xy_min.y / self.cell_size) as i32;

        let i_max = (xy_max.x / self.cell_size) as i32;
        let j_max = (xy_max.y / self.cell_size) as i32;

        for i in i_min..=i_max {
            for j in j_min..=j_max {
                self.cells
                    .entry((i, j))
                    .and_modify(|e| e.push(id))
                    .or_insert(vec![id]);
            }
        }

        self.objects.insert(id, obj);
    }

    fn remove(&mut self, id: T::Id) {
        if let Some(obj) = self.objects.remove(&id) {
            let Aabb { xy_min, xy_max } = obj.aabb();
            let i_min = (xy_min.x / self.cell_size) as i32;
            let j_min = (xy_min.y / self.cell_size) as i32;

            let i_max = (xy_max.x / self.cell_size) as i32;
            let j_max = (xy_max.y / self.cell_size) as i32;

            for i in i_min..=i_max {
                for j in j_min..=j_max {
                    if let Some(objs_in_cell) = self.cells.get_mut(&(i, j)) {
                        if let Some(index) = objs_in_cell.iter().position(|id_obj| *id_obj == id) {
                            objs_in_cell.swap_remove(index);
                        }
                    }
                }
            }
        }
    }

    fn update(&mut self, obj: T) {
        self.remove(obj.id());
        self.insert(obj);
    }

    fn query_point(&self, p: Point2<f32>) -> impl Iterator<Item = &T> {
        let i = (p.x / self.cell_size) as i32;
        let j = (p.y / self.cell_size) as i32;

        let objects = &self.objects;

        self.cells
            .get(&(i, j))
            .into_iter()
            .flat_map(move |ids| ids.iter().filter_map(move |id| objects.get(id)))
    }

    fn query_aabb(&self, region: Aabb) -> impl Iterator<Item = &T> {
        let Aabb { xy_min, xy_max } = region;
        let i_min = (xy_min.x / self.cell_size) as i32;
        let j_min = (xy_min.y / self.cell_size) as i32;

        let i_max = (xy_max.x / self.cell_size) as i32;
        let j_max = (xy_max.y / self.cell_size) as i32;

        let objects = &self.objects;
        let cells = &self.cells;

        (i_min..=i_max).flat_map(move |i| {
            (j_min..=j_max).flat_map(move |j| {
                cells
                    .get(&(i, j))
                    .into_iter()
                    .flat_map(move |ids| ids.iter().filter_map(move |id| objects.get(id)))
            })
        })
    }
}
