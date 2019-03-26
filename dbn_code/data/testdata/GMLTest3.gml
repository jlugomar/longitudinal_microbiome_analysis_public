<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
http://graphml.graphdrawing.org/xmlns/1.1/graphml.xsd">
<key id="key_weight" for="edge" attr.name="weight" attr.type="double"/>
<graph id="Full_exhaustive" edgedefault="directed">
<node id="class"/>
<node id="phenols"/>
<node id="flavanoids"/>
<node id="color"/>
<node id="hue"/>
<node id="od280od315"/>
<edge id="e1" source="class" target="phenols">
<data key="key_weight">59.225192</data>
</edge>
<edge id="e2" source="class" target="color">
<data key="key_weight">81.722275</data>
</edge>
<edge id="e3" source="class" target="hue">
<data key="key_weight">64.927116</data>
</edge>
<edge id="e4" source="class" target="od280od315">
<data key="key_weight">87.677879</data>
</edge>
<edge id="e5" source="phenols" target="flavanoids">
<data key="key_weight">104.338464</data>
</edge>
<edge id="e6" source="flavanoids" target="color">
<data key="key_weight">13.093560</data>
</edge>
<edge id="e7" source="hue" target="color">
<data key="key_weight">3.956640</data>
</edge>
<edge id="e8" source="od280od315" target="flavanoids">
<data key="key_weight">16.245673</data>
</edge>
</graph>
</graphml>
