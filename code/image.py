from math import floor
import os
import rawpy
import numpy as np
from PIL import Image, ImageEnhance
import cv2

def open_and_read_raw(file_path):
    try:
        with rawpy.imread(file_path) as raw:
            print(raw.raw_type)
            print(raw.color_desc)
            print(raw.raw_image.shape)
            print(raw.raw_image.dtype)

            raw_image = raw.raw_image.copy()  # Make a copy to avoid modifying original
            ob_value = determine_optical_black(raw)
            
           
    except Exception as e:
        print(f"Unable to open file: {file_path}. Error: {e}")
        return None, None

    return raw_image, ob_value


def determine_optical_black(raw):
    # Example: Extract optical black value from metadata
    ob_value = raw.black_level_per_channel  
    print("black_level_per_channel: ", ob_value)

    return ob_value
    
def optical_black(image, ob_value):
    # TODO: performance
    old_dtype = image.dtype
    print("raw bitdepth recognized as", old_dtype)
    res = image.astype(np.int64) - ob_value
    res[res < 0] = 0
    return res.astype(old_dtype)

def optical_black_per_channel(image, ob_value):
    # Subtract the optical black value from each channel
    image_corrected = np.zeros_like(image, dtype=np.int64)
    
    for c in range(3):
        image_corrected[:, :, c] = image[:, :, c].astype(np.int64) - ob_value[c]
        image_corrected[:, :, c][image_corrected[:, :, c] < 0] = 0
    
    return image_corrected.astype(image.dtype)

def white_balance(image):
    # TODO: probably wrong double check
    img_float = image.astype(np.float64)
    max_val = np.max(np.iinfo(image.dtype).max)
    print("max_val: ", max_val)
    avg_val = int(np.mean(img_float))
    print("avg_val: ", avg_val)
    scaled_img = img_float / (avg_val / np.iinfo(image.dtype).max)
    scaled_img = np.clip(scaled_img, np.iinfo(image.dtype).min, np.iinfo(image.dtype).max)
    
    return scaled_img.astype(image.dtype)

def white_balance_per_channel(image):
    img_float = image.astype(np.float64)
    max_val = np.iinfo(image.dtype).max
    
    # Compute scaling factors per channel based on the average value of each channel
    avg_val_per_channel = np.mean(img_float, axis=(0, 1))
    print("avg_val_per_channel: ", avg_val_per_channel, max_val)
    
    scaled_img = np.zeros_like(img_float)
    
    for c in range(3):
        scaling_factor = avg_val_per_channel[c] / max_val
        print("scaling_factor: ", scaling_factor)
        scaled_img[:, :, c] = img_float[:, :, c] / scaling_factor
        
        # Clip the values to the valid range
        scaled_img[:, :, c] = np.clip(scaled_img[:, :, c], np.iinfo(image.dtype).min, np.iinfo(image.dtype).max)
    
    return scaled_img.astype(image.dtype)


def demosaic_bayer_nearest(raw_image):
    """
    Convert a raw Bayer pattern image to an RGB image using nearest neighbor interpolation.
    
    Parameters:
    raw_image (np.ndarray): Input raw image of shape (height, width) with uint16 values.
    
    Returns:
    np.ndarray: RGB image of shape (height, width, 3) with uint16 values.
    """
    height, width = raw_image.shape
    rgb_image = np.zeros((height, width, 3), dtype=np.uint16)

    # Bayer pattern arrangement (assuming RGGB pattern)
    # R G
    # G B
    
    # Extracting the color channels
    red = raw_image[0::2, 0::2]
    green_red = raw_image[0::2, 1::2]
    green_blue = raw_image[1::2, 0::2]
    blue = raw_image[1::2, 1::2]

    # Initialize the channels
    red_channel = np.zeros((height, width), dtype=np.uint16)
    green_channel = np.zeros((height, width), dtype=np.uint16)
    blue_channel = np.zeros((height, width), dtype=np.uint16)

    # Assign the known values
    red_channel[0::2, 0::2] = red
    green_channel[0::2, 1::2] = green_red
    green_channel[1::2, 0::2] = green_blue
    blue_channel[1::2, 1::2] = blue

    # Nearest neighbor interpolation
    # For each position, replicate the nearest known value

    # Red channel
    red_channel[0::2, 1::2] = red_channel[0::2, 0::2]  # Right neighbors
    red_channel[1::2, 0::2] = red_channel[0::2, 0::2]  # Bottom neighbors
    red_channel[1::2, 1::2] = red_channel[0::2, 0::2]  # Bottom-right neighbors

    # Blue channel
    blue_channel[0::2, 1::2] = blue_channel[1::2, 1::2]  # Top neighbors
    blue_channel[1::2, 0::2] = blue_channel[1::2, 1::2]  # Left neighbors
    blue_channel[0::2, 0::2] = blue_channel[1::2, 1::2]  # Top-left neighbors

    # Green channel
    green_channel[0::2, 0::2] = green_channel[0::2, 1::2]  # Left neighbors
    green_channel[1::2, 1::2] = green_channel[1::2, 0::2]  # Right neighbors

    # Handle the edges for all channels
    red_channel[-1, :] = red_channel[-2, :]  # Last row
    red_channel[:, -1] = red_channel[:, -2]  # Last column

    blue_channel[-1, :] = blue_channel[-2, :]  # Last row
    blue_channel[:, -1] = blue_channel[:, -2]  # Last column

    green_channel[-1, :] = green_channel[-2, :]  # Last row
    green_channel[:, -1] = green_channel[:, -2]  # Last column

    # Combine channels into RGB image
    rgb_image[:, :, 0] = red_channel
    rgb_image[:, :, 1] = green_channel
    rgb_image[:, :, 2] = blue_channel

    return rgb_image


def calculate_gradients(image):
    grad_h = np.abs(image[:, 1:] - image[:, :-1])
    grad_v = np.abs(image[1:, :] - image[:-1, :])
    
    return grad_h, grad_v

def interpolate_channel(channel, mask, grad_h, grad_v):
    height, width = channel.shape
    interpolated = np.copy(channel)
    
    for y in range(height):
        for x in range(width):
            if not mask[y, x]:
                if x > 0 and x < width - 1:
                    left = channel[y, x - 1]
                    right = channel[y, x + 1]
                    grad_left = grad_h[y, x - 1]
                    grad_right = grad_h[y, x]
                    if grad_left < grad_right:
                        interpolated[y, x] = left
                    else:
                        interpolated[y, x] = right
                if y > 0 and y < height - 1:
                    top = channel[y - 1, x]
                    bottom = channel[y + 1, x]
                    grad_top = grad_v[y - 1, x]
                    grad_bottom = grad_v[y, x]
                    if grad_top < grad_bottom:
                        interpolated[y, x] = top
                    else:
                        interpolated[y, x] = bottom
    
    return interpolated

def demosaic_bayer_edge_filtered(raw_image):
    """
    Convert a raw Bayer pattern image to an RGB image using edge-filtered color difference interpolation.
    
    Parameters:
    raw_image (np.ndarray): Input raw image of shape (height, width) with uint16 values.
    
    Returns:
    np.ndarray: RGB image of shape (height, width, 3) with uint16 values.
    """
    height, width = raw_image.shape
    rgb_image = np.zeros((height, width, 3), dtype=np.uint16)

    # Bayer pattern arrangement (assuming RGGB pattern)
    # R G
    # G B

    # Extracting the color channels
    red = np.zeros((height, width), dtype=np.uint16)
    green = np.zeros((height, width), dtype=np.uint16)
    blue = np.zeros((height, width), dtype=np.uint16)

    red[0::2, 0::2] = raw_image[0::2, 0::2]
    green[0::2, 1::2] = raw_image[0::2, 1::2]
    green[1::2, 0::2] = raw_image[1::2, 0::2]
    blue[1::2, 1::2] = raw_image[1::2, 1::2]

    # Masks indicating known values
    mask_red = (red > 0)
    mask_green = (green > 0)
    mask_blue = (blue > 0)

    # Calculate gradients
    grad_h_red, grad_v_red = calculate_gradients(red)
    grad_h_green, grad_v_green = calculate_gradients(green)
    grad_h_blue, grad_v_blue = calculate_gradients(blue)

    # Interpolate missing values
    red = interpolate_channel(red, mask_red, grad_h_red, grad_v_red)
    green = interpolate_channel(green, mask_green, grad_h_green, grad_v_green)
    blue = interpolate_channel(blue, mask_blue, grad_h_blue, grad_v_blue)

    # Combine channels into RGB image
    rgb_image[:, :, 0] = red
    rgb_image[:, :, 1] = green
    rgb_image[:, :, 2] = blue

    return rgb_image



def apply_gamma_correction(image, gamma):
    # Normalize the image to the range [0, 1]
    image_normalized = image / 65535.0
    
    # Apply gamma correction
    image_gamma_corrected = np.power(image_normalized, gamma)
    
    # Scale back to the range [0, 65535]
    image_corrected = np.uint16(image_gamma_corrected * 65535)
    
    return image_corrected


def show_images(raw_image, processed_image):
    raw_bgr = cv2.cvtColor(raw_image, cv2.COLOR_RGB2BGR)
    processed_bgr = cv2.cvtColor(processed_image, cv2.COLOR_RGB2BGR)
    
    # Resize images to fit Full HD (1920x1080)
    raw_bgr_resized = cv2.resize(raw_bgr, (1920, 1080))
    processed_bgr_resized = cv2.resize(processed_bgr, (1920, 1080))
    
    # raw_bgr_resized = raw_bgr
    # processed_bgr_resized = processed_bgr

    current_image = raw_bgr_resized
    cv2.imshow('Image Comparison', current_image)
    
    while True:
        key = cv2.waitKey(1) & 0xFF
        if key == ord(' '):  # Space key pressed
            current_image = processed_bgr_resized if (current_image is raw_bgr_resized) else raw_bgr_resized
            cv2.imshow('Image Comparison', current_image)
        elif key == 27:  # Escape key pressed
            break
    
    cv2.destroyAllWindows()

if __name__ == '__main__':
    file_path = r"/home/aryan/DigitalCameraContest/1/SRF/DSC00474.SRF"
    output_path = r"/home/aryan/DigitalCameraContest/1/process.jpg"

    # Open and read raw image data
    raw_image, ob_value = open_and_read_raw(file_path)
    
    if raw_image is not None and ob_value is not None:
        processed_image = raw_image

        processed_image = optical_black(raw_image, ob_value[0])
        processed_image = white_balance(processed_image)
        # TODO: debatable whether to apply color interpolation before or after white balance
        # if wb first then gray world on grayscale, otherwise per channel coefficients see awb.ppt
        
        processed_image = demosaic_bayer_nearest(processed_image)
        # processed_image = demosaic_bayer_edge_filtered(processed_image)

        # processed_image = white_balance_per_channel(processed_image)
        # processed_image = optical_black_per_channel(processed_image, ob_value[:3])
        processed_image = apply_gamma_correction(processed_image, 0.45)


        # processed_image = color_correction(processed_image)
        # processed_image = gamma_correction(processed_image)

        # processed_image = convert_color_space(processed_image, 'sRGB1', 'CAM02-UCS')
        # processed_image = adjust_hue_saturation(processed_image, hue=0, saturation=1.2)
        # processed_image = convert_color_space(processed_image, 'CAM02-UCS', 'sRGB1')
        # processed_image = edge_enhancement(processed_image)

        show_images(raw_image, processed_image)
        
        # Save the processed image (optional)
        # image_compression(processed_image, output_path)
        print(f"Processed image saved to: {output_path}")