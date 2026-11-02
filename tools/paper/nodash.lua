-- The author's rule: no en or em dashes anywhere. citeproc writes page ranges (from "--" in the .bib) and runs of
-- citations (superscript 1 to 3) with an en dash; the first pass turns those into hyphens. The second pass does the
-- same for any dash left elsewhere (for example "--" in the Markdown, which pandoc's smart typography turns into an en dash)
-- and warns, since a dash in the prose should be rewritten, not patched.
local function fix(s) return (s:gsub("\u{2013}", "-"):gsub("\u{2014}", "-")) end

local quiet = { Str = function(s) return pandoc.Str(fix(s.text)) end }

-- Pass 1: the reference list, and citations (Nature style writes runs of citations such as 1 to 3 with an en dash).
local refs = {
  Div = function(d)
    if d.identifier == "refs" then return pandoc.walk_block(d, quiet) end
  end,
  Cite = function(c)
    c.content = pandoc.walk_inline(pandoc.Span(c.content), quiet).content
    return c
  end,
}

-- Pass 2: anything else.
local prose = {
  Str = function(s)
    local t = fix(s.text)
    if t ~= s.text then
      io.stderr:write("nodash.lua: dash outside the references, rewrite it: " .. s.text .. "\n")
      return pandoc.Str(t)
    end
  end,
}

return { refs, prose }
