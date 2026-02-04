#!/usr/bin/env python3
"""Script to update color from dark purple to lighter indigo"""

with open('app.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Replace dark purple with lighter indigo
old_color = '#667eea'
new_color = '#6366f1'
count = content.count(old_color)

content = content.replace(old_color, new_color)

with open('app.py', 'w', encoding='utf-8') as f:
    f.write(content)

print(f"✅ Successfully replaced {count} instances of {old_color} with {new_color}")
print(f"🎨 New color: Lighter Indigo ({new_color}) - brighter and less heavy!")
