```@eval
using Markdown
# Read the README.md file as a string
readme = read("../../README.md", String)

# Find the position of the first '# GeMotion.jl' character (first headline)
pos = findfirst(r"# GeMotion.jl", readme)[1]

# Check if a '# GeMotion.jl' was found
if isnothing(pos)
    error("No '# GeMotion.jl' found in README.md")
end

print(readme)
print(pos)

# Extract the content from the first headline onward
readme_content = readme[pos:end]

# Parse the extracted content as Markdown and include it
Markdown.parse(readme_content)
```

## API

- See [API reference](./api.md)
