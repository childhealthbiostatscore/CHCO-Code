function Image(img)
  -- Return the image with the path appended as part of the caption
  local caption = img.caption
  table.insert(caption, pandoc.LineBreak())
  table.insert(caption, pandoc.Str("Path: "))
  table.insert(caption, pandoc.Code(img.src))
  
  return pandoc.Image(caption, img.src, img.title, img.attr)
end