#![crate_type = "lib"]
#![cfg_attr(not(test), no_std)]

use ringbuffer::{ConstGenericRingBuffer, RingBuffer, RingBufferExt, RingBufferWrite};

use embedded_graphics::prelude::*;
use embedded_graphics::primitives::{Line, PrimitiveStyle, Rectangle};

use time_series_filter::{EwmaFilter, FloatSeriesEwmaFilter };


pub struct Emplot<C, const N: usize>
    where
        C: PixelColor,
{
    /// stores N samples
    samples: ConstGenericRingBuffer<f32, N>,
    /// outer pixel bounds of this plot
    bounds: Rectangle,
    /// defines the number of most recent values that we'll draw
    n_draw_samples: usize,
    /// color to use for drawing the plot itself
    line_color: C,
    /// width of the stroke to use for drawing the plot iself
    line_stroke_w: u32,

    /// cached maximum pixel height for the plot line
    max_pixel_height: u32,
 
    iir_filter: FloatSeriesEwmaFilter<f32>,

}

impl<C, const N: usize> Emplot<C, N>
    where
        C: PixelColor,
{
    pub fn new(
        bounds: Rectangle,
        n_draw_samples: usize,
        line_color: C,
        line_stroke_width: u32,
    ) -> Self {
        assert!(N >= n_draw_samples, "n_draw_samples can't exceed capacity");
        Self {
            samples: ConstGenericRingBuffer::<f32, N>::new(),
            n_draw_samples,
            bounds,
            line_color,
            line_stroke_w: line_stroke_width,
            max_pixel_height: (bounds.top_left.y + bounds.size.height as i32) as u32,
            iir_filter: FloatSeriesEwmaFilter::default(),
        }
    }

    /// Push a new sample value to the sample buffer
    pub fn push(&mut self, val: f32) {
        self.samples.push(val);
        self.iir_filter.push_sample(val);
    }

    /// The number of samples we will attempt to draw, currently
    pub fn num_draw_items(&self) -> usize {
        let n_items =
            if self.samples.len() > self.n_draw_samples { self.n_draw_samples}
            else { self.samples.len() };
        n_items
    }

}

impl<C, const N: usize> Drawable for Emplot<C, N>
    where
        C: PixelColor,
{
    type Color = C;
    type Output = ();

    fn draw<D>(&self, target: &mut D) -> Result<Self::Output, D::Error>
        where
            D: DrawTarget<Color=Self::Color>,
    {
        let plot_style = PrimitiveStyle::with_stroke(self.line_color, self.line_stroke_w);
        let mut slope: f32 = (self.bounds.size.height - self.line_stroke_w) as f32;
        // the sample index at which we should start drawing
        // may be different from the first sample
        let start_idx =
            if self.samples.len() > self.n_draw_samples {
                self.samples.len() - self.n_draw_samples
            }
            else {0};

        /*
        // find the minimum and maximum samples in the entire sample buffer
        let (min, max): (&f32, &f32) =
            self.samples.iter()
                .fold((&f32::MAX, &f32::MIN), |mut acc, val| {
                    if val < acc.0 { acc.0 = val; }
                    if val > acc.1 { acc.1 = val; }
                    acc
                });
        */

        // use exponential weighted moving average for recent min/max
        let range = self.iir_filter.local_range();
        let display_min = range.start;
	let display_max = range.end;


        // recalculate slope based on min and max
        if display_max != display_min {
            slope /= display_max - display_min;
        }

        // the samples len may be shorter than n_draw_samples, initially
        let n_items = self.num_draw_items();
        if n_items < 1 { return Ok(()) };
        let px_per_seg = (self.bounds.size.width - 1) as f32 / (n_items - 1) as f32;
        let mut last_pt = Point::new(0, 0);

        for (i, val) in self.samples.iter().skip(start_idx).enumerate() {
            let idx_offset = if i > start_idx { i - start_idx } else { 0 };
            let scaled_val =
                self.max_pixel_height as f32
                - ((val - display_min) as f32 * slope)
                - self.line_stroke_w as f32 / 2f32;

            let pt_x = (idx_offset as f32 * px_per_seg) as i32 + self.bounds.top_left.x;
            let pt_y = scaled_val as i32;
            let cur_pt = Point::new(pt_x, pt_y);

            // skip first point as it goes from zero
            if idx_offset > 0  {
                // draw using supplied closure drawing function
                Line::new(last_pt, cur_pt)
                    .into_styled( plot_style)
                    .draw(target)?;
            }
            last_pt = cur_pt;
        }

        Ok(())

    }
}


#[cfg(test)]
mod tests {
    use embedded_graphics::pixelcolor::BinaryColor;
    use embedded_graphics::primitives::{Rectangle};
    use super::*;

    #[test]
    fn push_data() {
        const MAX_SAMPLES: usize = 64;
        const MAX_DRAW_ITEMS: usize = 32;
        let bounds = Rectangle::new(Point::new(0, 0), Size::new(120, 120));
        let mut emplot = Emplot::<_,MAX_SAMPLES>::new(
            bounds,
            MAX_DRAW_ITEMS,
            BinaryColor::On,
            1,
        );

        for i in 0..3 { emplot.push(i as f32);}
        assert_eq!(emplot.num_draw_items(), 3);

        for i in 3..(2*MAX_DRAW_ITEMS) { emplot.push(i as f32);}
        assert_eq!(emplot.num_draw_items(), MAX_DRAW_ITEMS);

    }
}
