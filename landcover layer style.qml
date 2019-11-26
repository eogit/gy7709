<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis styleCategories="AllStyleCategories" version="3.8.2-Zanzibar" hasScaleBasedVisibilityFlag="0" minScale="1e+08" maxScale="0">
  <flags>
    <Identifiable>1</Identifiable>
    <Removable>1</Removable>
    <Searchable>1</Searchable>
  </flags>
  <customproperties>
    <property value="false" key="WMSBackgroundLayer"/>
    <property value="false" key="WMSPublishDataSourceUrl"/>
    <property value="0" key="embeddedWidgets/count"/>
    <property value="Value" key="identify/format"/>
  </customproperties>
  <pipe>
    <rasterrenderer band="1" alphaBand="-1" opacity="1" type="paletted">
      <rasterTransparency/>
      <minMaxOrigin>
        <limits>None</limits>
        <extent>WholeRaster</extent>
        <statAccuracy>Estimated</statAccuracy>
        <cumulativeCutLower>0.02</cumulativeCutLower>
        <cumulativeCutUpper>0.98</cumulativeCutUpper>
        <stdDevFactor>2</stdDevFactor>
      </minMaxOrigin>
      <colorPalette>
        <paletteEntry alpha="255" value="0" color="#000004" label="Missing"/>
        <paletteEntry alpha="255" value="1" color="#1face7" label="Water"/>
        <paletteEntry alpha="255" value="2" color="#bd282a" label="Residential"/>
        <paletteEntry alpha="255" value="3" color="#a30e72" label="Industrial"/>
        <paletteEntry alpha="255" value="4" color="#b8dba3" label="Pasture"/>
        <paletteEntry alpha="255" value="5" color="#ffd09a" label="Crops"/>
        <paletteEntry alpha="255" value="6" color="#c68567" label="Bare soil"/>
        <paletteEntry alpha="255" value="7" color="#728454" label="Forest"/>
      </colorPalette>
      <colorramp name="[source]" type="random">
        <prop k="count" v="10"/>
        <prop k="hueMax" v="359"/>
        <prop k="hueMin" v="0"/>
        <prop k="rampType" v="random"/>
        <prop k="satMax" v="240"/>
        <prop k="satMin" v="100"/>
        <prop k="valMax" v="240"/>
        <prop k="valMin" v="200"/>
      </colorramp>
    </rasterrenderer>
    <brightnesscontrast contrast="0" brightness="0"/>
    <huesaturation colorizeRed="255" colorizeStrength="100" colorizeGreen="128" colorizeOn="0" colorizeBlue="128" grayscaleMode="0" saturation="0"/>
    <rasterresampler maxOversampling="2"/>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
