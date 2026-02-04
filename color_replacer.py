#!/usr/bin/env python3
import re

# Read the file
with open('app.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Count occurrences
old_color = '#667eea'
new_color = '#6366f1'
count_before = content.count(old_color)

# Replace all instances
content = content.replace(old_color, new_color)

# Also update the rgba references
content = content.replace('rgba(102, 126, 234', 'rgba(99, 102, 241')

# Write back
with open('app.py', 'w', encoding='utf-8') as f:
    f.write(content)

count_after = content.count(old_color)

print(f"✅ Color replacement complete!")
print(f"📊 Replaced {count_before - count_after} instances")
print(f"🎨 Old color: {old_color} (dark purple - heavy)")  
print(f"🎨 New color: {new_color} (lighter indigo - bright and modern!)")
print(f"✨ Remaining instances of old color: {count_after}")
