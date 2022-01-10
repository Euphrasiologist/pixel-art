use draw_prim::{self};
use image::{Rgb, RgbImage};
fn main() {
    let mut canvas = RgbImage::new(200, 200);

    for (_x, _y, pixel) in canvas.enumerate_pixels_mut() {
        let curr_color = Rgb([223, 30, 30]); // hardcoded background color
        *pixel = curr_color;
    }

    // draw_prim::draw_line_segment_mut(&mut canvas, (0.0, 0.0), (54.0, 90.0), Rgb([10, 100, 100]));
    draw_prim::draw_filled_circle_mut(&mut canvas, (80, 30), 25, Rgb([245, 255, 255]));

    let telegraph_data = [(10.0, 50.0), (50.0, 50.0), (70.0, 90.0), (170.0, 90.0)].windows(2);

    for el in telegraph_data {
        println!("{:?} {:?}", el[0], el[1]);
        draw_prim::draw_catenary_curve(&mut canvas, el[0], el[1], Rgb([0, 0, 0]), 0.0);
    }

    let rect = draw_prim::Rect::at(0, 0).of_size(200, 200);

    draw_prim::draw_filled_rect_mut(&mut canvas, rect, Rgb([250, 250, 250]));

    canvas.save("output.png").unwrap();
}
