import os
import rawpy
import numpy as np
from PIL import Image, ImageEnhance
import cv2

def open_and_read_raw(file_path):
    try:
        with rawpy.imread(file_path) as raw:
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


def white_balance(image):
    img_float = image.astype(np.float64)
    max_val = np.max(img_float)
    print("max_val: ", max_val)
    avg_val = int(np.mean(img_float))
    print("avg_val: ", avg_val)
    scaled_img = img_float / (avg_val / np.iinfo(image.dtype).max)
    scaled_img = np.clip(scaled_img, np.iinfo(image.dtype).min, np.iinfo(image.dtype).max)
    
    return scaled_img.astype(image.dtype)

def color_interpolation(image):
    # Placeholder for color interpolation
    return image

def color_correction(image):
    # Placeholder for color correction
    return image

def gamma_correction(image):
    # Placeholder for gamma correction
    return image

def convert_color_space(image, from_space, to_space):
    # Placeholder for color space conversion
    return image

def adjust_hue_saturation(image, hue=0, saturation=1.0):
    # Placeholder for hue/saturation adjustment
    return image

def edge_enhancement(rgb_image, space='RGB', strength=1.0):
    if space == 'YUV':
        
        yuv_image = cv2.cvtColor(rgb_image, cv2.COLOR_RGB2YUV)
        y, u, v = cv2.split(yuv_image)

        y_enhanced = cv2.equalizeHist(y)
        y_enhanced = cv2.GaussianBlur(y_enhanced, (5, 5), 0)
        y_enhanced = cv2.addWeighted(y, 1.0 + strength, y_enhanced, -strength, 0)
        
        enhanced_yuv_image = cv2.merge([y_enhanced, u, v])

        enhanced_rgb_image = cv2.cvtColor(enhanced_yuv_image, cv2.COLOR_YUV2RGB)
        
        return enhanced_rgb_image
    
    elif space == 'RGB': 
        image = rgb_image.astype(np.float32)
        
        
        kernel = np.array([[-1.0, -1, -1],
                       [-1,  9, -1],
                       [-1, -1, -1]])
        
        kernel *= strength
        
        enhanced_image = np.zeros_like(image)
        for channel in range(image.shape[2]):
            enhanced_image[:,:,channel] = ndimage.convolve(image[:,:,channel], kernel)

        enhanced_image = np.clip(enhanced_image, 0, 255)

        enhanced_image = enhanced_image.astype(np.uint8)
        
        return enhanced_image
    else:
        print ('Not a valid space, skipping')     
        return enhanced_image

def image_compression(image, output_path):
    compressed_image = Image.fromarray(image)
    compressed_image.save(output_path, 'JPEG', quality=90)

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
    file_path = r"C:\Users\IDF\OneDrive\taw\code\subject\1.SRF"
    output_path = r"C:\Users\IDF\OneDrive\taw\code\subject\processed_image.jpg"

    # Open and read raw image data
    raw_image, ob_value = open_and_read_raw(file_path)
    
    if raw_image is not None and ob_value is not None:
        processed_image = raw_image

        processed_image = optical_black(raw_image, ob_value[0])

        processed_image = white_balance(processed_image)
        # processed_image = color_interpolation(processed_image)
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