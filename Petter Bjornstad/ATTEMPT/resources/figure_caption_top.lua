-- Counter for figures
local fig_counter = 0

function Image(img)
  -- Only process images that have captions
  if img.caption and #img.caption > 0 then
    fig_counter = fig_counter + 1
    
    -- Get caption text
    local caption_text = pandoc.utils.stringify(img.caption)
    
    -- Create new caption with "Figure N: " prefix
    local new_caption = {}
    table.insert(new_caption, pandoc.Strong{pandoc.Str("Figure " .. fig_counter .. ": ")})
    
    -- Add the rest of the caption
    for _, elem in ipairs(img.caption) do
      table.insert(new_caption, elem)
    end
    
    -- Update the image caption
    img.caption = new_caption
    
    return img
  end
end

-- For images inside Para elements
function Para(para)
  local has_image = false
  local new_content = {}
  
  for i, elem in ipairs(para.content) do
    if elem.t == "Image" and elem.caption and #elem.caption > 0 then
      has_image = true
      fig_counter = fig_counter + 1
      
      -- Create caption paragraph
      local caption_para = pandoc.Para{
        pandoc.Strong{pandoc.Str("Figure " .. fig_counter .. ": " .. pandoc.utils.stringify(elem.caption))}
      }
      
      -- Clear the image caption
      elem.caption = {}
      
      -- Return both caption and image
      return {caption_para, pandoc.Para{elem}}
    end
  end
  
  return para
end