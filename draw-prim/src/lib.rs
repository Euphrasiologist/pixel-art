// we need a canvas
// and a bunch of drawing primitives
// lines, paths, arcs,
// circles, rectangles, ellipses

// almost all poached from imageproc
// https://docs.rs/imageproc/

use image::{GenericImage, GenericImageView, Pixel, Rgb};
use std::mem::swap;

pub trait Canvas {
    /// The type of `Pixel` that can be drawn on this canvas.
    type Pixel: Pixel;

    /// The width and height of this canvas.
    fn dimensions(&self) -> (u32, u32);

    /// The width of this canvas.
    fn width(&self) -> u32 {
        self.dimensions().0
    }

    /// The height of this canvas.
    fn height(&self) -> u32 {
        self.dimensions().1
    }

    /// Returns the pixel located at (x, y).
    fn get_pixel(&self, x: u32, y: u32) -> Self::Pixel;

    /// Draw a pixel at the given coordinates. `x` and `y`
    /// should be within `dimensions` - if not then panicking
    /// is a valid implementation behaviour.
    fn draw_pixel(&mut self, x: u32, y: u32, colour: Self::Pixel);
}

impl<I> Canvas for I
where
    I: GenericImage,
{
    type Pixel = I::Pixel;

    fn dimensions(&self) -> (u32, u32) {
        <I as GenericImageView>::dimensions(self)
    }

    fn get_pixel(&self, x: u32, y: u32) -> Self::Pixel {
        self.get_pixel(x, y)
    }

    fn draw_pixel(&mut self, x: u32, y: u32, colour: Self::Pixel) {
        self.put_pixel(x, y, colour)
    }
}

// primitives

/// Draws as much of the line segment between start and end as lies inside the image bounds.
/// Uses [Bresenham's line drawing algorithm](https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm).
pub fn draw_line_segment_mut<C>(
    canvas: &mut C,
    start: (f32, f32),
    end: (f32, f32),
    colour: C::Pixel,
) where
    C: Canvas,
    C::Pixel: 'static,
{
    let (width, height) = canvas.dimensions();
    let in_bounds = |x, y| x >= 0 && x < width as i32 && y >= 0 && y < height as i32;

    let line_iterator = BresenhamLineIter::new(start, end);

    for point in line_iterator {
        let x = point.0;
        let y = point.1;

        if in_bounds(x, y) {
            canvas.draw_pixel(x as u32, y as u32, colour);
        }
    }
}

pub struct BresenhamLineIter {
    dx: f32,
    dy: f32,
    x: i32,
    y: i32,
    error: f32,
    end_x: i32,
    is_steep: bool,
    y_step: i32,
}

impl BresenhamLineIter {
    /// Creates a [`BresenhamLineIter`](struct.BresenhamLineIter.html) which will iterate over the integer coordinates
    /// between `start` and `end`.
    pub fn new(start: (f32, f32), end: (f32, f32)) -> BresenhamLineIter {
        let (mut x0, mut y0) = (start.0, start.1);
        let (mut x1, mut y1) = (end.0, end.1);

        let is_steep = (y1 - y0).abs() > (x1 - x0).abs();
        if is_steep {
            swap(&mut x0, &mut y0);
            swap(&mut x1, &mut y1);
        }

        if x0 > x1 {
            swap(&mut x0, &mut x1);
            swap(&mut y0, &mut y1);
        }

        let dx = x1 - x0;

        BresenhamLineIter {
            dx,
            dy: (y1 - y0).abs(),
            x: x0 as i32,
            y: y0 as i32,
            error: dx / 2f32,
            end_x: x1 as i32,
            is_steep,
            y_step: if y0 < y1 { 1 } else { -1 },
        }
    }
}

impl Iterator for BresenhamLineIter {
    type Item = (i32, i32);

    fn next(&mut self) -> Option<(i32, i32)> {
        if self.x > self.end_x {
            None
        } else {
            let ret = if self.is_steep {
                (self.y, self.x)
            } else {
                (self.x, self.y)
            };

            self.x += 1;
            self.error -= self.dy;
            if self.error < 0f32 {
                self.y += self.y_step;
                self.error += self.dx;
            }

            Some(ret)
        }
    }
}

// Draw as much of a circle, including its contents, as lies inside the image bounds.
pub fn draw_filled_circle_mut<C>(canvas: &mut C, center: (i32, i32), radius: i32, colour: C::Pixel)
where
    C: Canvas,
    C::Pixel: 'static,
{
    let mut x = 0i32;
    let mut y = radius;
    let mut p = 1 - radius;
    let x0 = center.0;
    let y0 = center.1;

    while x <= y {
        draw_line_segment_mut(
            canvas,
            ((x0 - x) as f32, (y0 + y) as f32),
            ((x0 + x) as f32, (y0 + y) as f32),
            colour,
        );
        draw_line_segment_mut(
            canvas,
            ((x0 - y) as f32, (y0 + x) as f32),
            ((x0 + y) as f32, (y0 + x) as f32),
            colour,
        );
        draw_line_segment_mut(
            canvas,
            ((x0 - x) as f32, (y0 - y) as f32),
            ((x0 + x) as f32, (y0 - y) as f32),
            colour,
        );
        draw_line_segment_mut(
            canvas,
            ((x0 - y) as f32, (y0 - x) as f32),
            ((x0 + y) as f32, (y0 - x) as f32),
            colour,
        );

        x += 1;
        if p < 0 {
            p += 2 * x + 1;
        } else {
            y -= 1;
            p += 2 * (x - y) + 1;
        }
    }
}

// from this JS implementation
// https://github.com/dulnan/catenary-curve/blob/5894b52b2c661a95218344cf5f4d02e61c360909/src/Catenary.js
pub fn draw_catenary_curve<C>(
    canvas: &mut C,
    start: (f32, f32),
    end: (f32, f32),
    colour: C::Pixel,
    mut chain_length: f32,
) where
    C: Canvas,
    C::Pixel: 'static,
{
    fn get_catenary_parameter(h: f32, v: f32, length: f32, limit: i32, EPSILON: f32) -> f32 {
        let m = (length * length - v * v).sqrt() / h;
        let mut x = (m).acosh() + 1.0;
        let mut prevx = -1.0;
        let mut count = 0;

        while (x - prevx).abs() > EPSILON && count < limit {
            prevx = x;
            x = x - ((x).sinh() - m * x) / ((x).cosh() - m);
            count += 1;
        }

        h / (2.0 * x)
    }

    fn get_curve(
        a: f32,
        p1: (f32, f32),
        p2: (f32, f32),
        offset_x: f32,
        offset_y: f32,
        segments: i32,
    ) -> Vec<(f32, f32)> {
        let mut data = Vec::new();
        data.push((p1.0, a * ((p1.0 - offset_x) / a).cosh() + offset_y));

        let d = p2.0 - p1.0;
        let length = segments - 1;

        for i in 0..length {
            let x = p1.0 + d * (i as f32 + 0.5) / length as f32;
            let y = a * ((x - offset_x) / a).cosh() + offset_y;
            data.push((x, y));
        }

        data.push((p2.0, a * ((p2.0 - offset_x) / a).cosh() + offset_y));

        data
    }

    let EPSILON = 1e-6;
    let segments = 50;
    let iteration_limit = 100;

    // so we can flip
    let p1_ = start;
    let p2_ = end;

    let is_flipped = start.0 > end.0;

    let p1: (f32, f32);
    let p2: (f32, f32);

    match is_flipped {
        true => {
            p1 = p2_;
            p2 = p1_;
        }
        false => {
            p1 = p1_;
            p2 = p2_;
        }
    }
    // get distance between points
    let diff = (p1.0 - p2.0, p1.1 - p2.1);
    let dist = (diff.0.powi(2) + diff.1.powi(2)).sqrt();

    if chain_length < dist {
        chain_length = dist + 10.0;
    }

    let curve_data;
    // let mut is_straight = true;

    // if dist < chain_length {
    let diff_x = p2.0 - p1.0;

    if diff_x > 0.01 {
        let h = p2.0 - p1.0;
        let v = p2.1 - p1.1;
        let a = -get_catenary_parameter(h, v, chain_length, iteration_limit, EPSILON);
        let x = (a * ((chain_length + v) / (chain_length - v)).log(10.0) - h) * 0.5;
        let y = a * (x / a).cosh();
        let offset_x = p1.0 - x;
        let offset_y = p1.1 - y;
        curve_data = get_curve(a, p1, p2, offset_x, offset_y, segments);
        // is_straight = false;
    } else {
        let mx = (p1.0 + p2.0) * 0.5;
        let my = (p1.1 + p2.1 + chain_length) * 0.5;

        curve_data = vec![(p1.0, p1.1), (mx, my), (p2.0, p2.1)];
    }
    // } else {
    //     curve_data = vec![(p1.0, p1.1), (p2.0, p2.1)];
    // }

    // now the data has been calculated
    let mut index = 0;
    loop {
        if index == curve_data.len() - 1 {
            break;
        }
        let t_start = curve_data[index];
        let t_end = curve_data[index + 1];

        // println!("start: {:?}, end: {:?}", t_start, t_end);

        draw_line_segment_mut(canvas, t_start, t_end, colour);

        index += 1;
    }
}

// Draws as much of a cubic bezier curve as lies within image bounds.
pub fn draw_cubic_bezier_curve_mut<C>(
    canvas: &mut C,
    start: (f32, f32),
    end: (f32, f32),
    control_a: (f32, f32),
    control_b: (f32, f32),
    colour: C::Pixel,
) where
    C: Canvas,
    C::Pixel: 'static,
{
    // Bezier Curve function from: https://pomax.github.io/bezierinfo/#control
    let cubic_bezier_curve = |t: f32| {
        let t2 = t * t;
        let t3 = t2 * t;
        let mt = 1.0 - t;
        let mt2 = mt * mt;
        let mt3 = mt2 * mt;
        let x = (start.0 * mt3)
            + (3.0 * control_a.0 * mt2 * t)
            + (3.0 * control_b.0 * mt * t2)
            + (end.0 * t3);
        let y = (start.1 * mt3)
            + (3.0 * control_a.1 * mt2 * t)
            + (3.0 * control_b.1 * mt * t2)
            + (end.1 * t3);
        (x.round(), y.round()) // round to nearest pixel, to avoid ugly line artifacts
    };

    let distance = |point_a: (f32, f32), point_b: (f32, f32)| {
        ((point_a.0 - point_b.0).powi(2) + (point_a.1 - point_b.1).powi(2)).sqrt()
    };

    // Approximate curve's length by adding distance between control points.
    let curve_length_bound: f32 =
        distance(start, control_a) + distance(control_a, control_b) + distance(control_b, end);

    // Use hyperbola function to give shorter curves a bias in number of line segments.
    let num_segments: i32 = ((curve_length_bound.powi(2) + 800.0).sqrt() / 8.0) as i32;

    // Sample points along the curve and connect them with line segments.
    let t_interval = 1f32 / (num_segments as f32);
    let mut t1 = 0f32;
    for i in 0..num_segments {
        let t2 = (i as f32 + 1.0) * t_interval;
        draw_line_segment_mut(
            canvas,
            cubic_bezier_curve(t1),
            cubic_bezier_curve(t2),
            colour,
        );
        t1 = t2;
    }
}

// Basic manipulation of rectangles.
use rand::prelude::*;

// Draw as much of a rectangle, including its boundary, as lies inside the image bounds.
pub fn draw_filled_rect_mut<C>(canvas: &mut C, rect: Rect, color: C::Pixel)
where
    C: Canvas,
    C::Pixel: 'static,
{
    let mut rng = thread_rng();

    let canvas_bounds = Rect::at(0, 0).of_size(canvas.width(), canvas.height());
    if let Some(intersection) = canvas_bounds.intersect(rect) {
        for dy in (0..intersection.height()).choose_multiple(&mut rng, 100) {
            for dx in (0..intersection.width()).choose_multiple(&mut rng, 5) {
                let x = intersection.left() as u32 + dx;
                let y = intersection.top() as u32 + dy;
                canvas.draw_pixel(x, y, color);
            }
        }
    }
}

use std::cmp;

/// A rectangular region of non-zero width and height.

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct Rect {
    left: i32,
    top: i32,
    width: u32,
    height: u32,
}

/// A geometrical representation of a set of 2D points with coordinate type T.
pub trait Region<T> {
    /// Whether this region contains the given point.
    fn contains(&self, x: T, y: T) -> bool;
}

impl Rect {
    /// Reduces possibility of confusing coordinates and dimensions
    /// when specifying rects.
    ///
    /// See the [struct-level documentation](struct.Rect.html) for examples.
    pub fn at(x: i32, y: i32) -> RectPosition {
        RectPosition { left: x, top: y }
    }

    /// Smallest y-coordinate reached by rect.
    ///
    /// See the [struct-level documentation](struct.Rect.html) for examples.
    pub fn top(&self) -> i32 {
        self.top
    }

    /// Smallest x-coordinate reached by rect.
    ///
    /// See the [struct-level documentation](struct.Rect.html) for examples.
    pub fn left(&self) -> i32 {
        self.left
    }

    /// Greatest y-coordinate reached by rect.
    ///
    /// See the [struct-level documentation](struct.Rect.html) for examples.
    pub fn bottom(&self) -> i32 {
        self.top + (self.height as i32) - 1
    }

    /// Greatest x-coordinate reached by rect.
    ///
    /// See the [struct-level documentation](struct.Rect.html) for examples.
    pub fn right(&self) -> i32 {
        self.left + (self.width as i32) - 1
    }

    /// Width of rect.
    pub fn width(&self) -> u32 {
        self.width
    }

    /// Height of rect.
    pub fn height(&self) -> u32 {
        self.height
    }

    /// Returns the intersection of self and other, or none if they are are disjoint.

    pub fn intersect(&self, other: Rect) -> Option<Rect> {
        let left = cmp::max(self.left, other.left);
        let top = cmp::max(self.top, other.top);
        let right = cmp::min(self.right(), other.right());
        let bottom = cmp::min(self.bottom(), other.bottom());

        if right < left || bottom < top {
            return None;
        }

        Some(Rect {
            left,
            top,
            width: (right - left) as u32 + 1,
            height: (bottom - top) as u32 + 1,
        })
    }
}

impl Region<i32> for Rect {
    fn contains(&self, x: i32, y: i32) -> bool {
        self.left <= x && x <= self.right() && self.top <= y && y <= self.bottom()
    }
}

impl Region<f32> for Rect {
    fn contains(&self, x: f32, y: f32) -> bool {
        self.left as f32 <= x
            && x <= self.right() as f32
            && self.top as f32 <= y
            && y <= self.bottom() as f32
    }
}

/// Position of the top left of a rectangle.
/// Only used when building a [`Rect`](struct.Rect.html).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct RectPosition {
    left: i32,
    top: i32,
}

impl RectPosition {
    /// Construct a rectangle from a position and size. Width and height
    /// are required to be strictly positive.
    ///
    /// See the [`Rect`](struct.Rect.html) documentation for examples.
    pub fn of_size(self, width: u32, height: u32) -> Rect {
        assert!(width > 0, "width must be strictly positive");
        assert!(height > 0, "height must be strictly positive");
        Rect {
            left: self.left,
            top: self.top,
            width,
            height,
        }
    }
}
