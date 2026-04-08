# PatchedConicFixes

This is a near overhaul of KSP 1.12.x's encounter resolution in the patched conic
solver.  It fixes several regressions introduced around KSP 1.8.0, along with
multiple design problems that KSP had from the start.  Many or most of the bugs
where SOI encounters to go undetected, even when your ship is on a clear intercept
course, have been fixed.  A few examples of affected scenarios include problems
with flickering orbits around Jool, and transfers to an ascending Moon in
RealSolarSystem, and many others. The replacement encounter solver is also
approximately 5× faster than the stock implementation.

## Requirements

- Kerbal Space Program 1.12.x

## Installation

Install via CKAN, or manually copy `GameData/PatchedConicFixes` into your KSP `GameData` folder.

## Building from Source

### Prerequisites

- [.NET SDK 6.0+](https://dotnet.microsoft.com/download)
- KSP 1.12.x installed

### Setup

Edit `Source/Directory.Build.props` to point to your KSP installation:

```xml
<KSPData>/path/to/KSP.app/Contents/Resources/Data</KSPData>     <!-- Mac -->
<KSPGameData>/path/to/KSP/GameData</KSPGameData>
```

Or create a `Source/Directory.Build.props.user` (gitignored) to keep your local paths out of version control:

```xml
<Project>
  <PropertyGroup>
    <KSPData>/your/local/path/KSP.app/Contents/Resources/Data</KSPData>
    <KSPGameData>/your/local/path/GameData</KSPGameData>
  </PropertyGroup>
</Project>
```

### Build

```
cd Source
dotnet build
```

The compiled DLL will be placed in `GameData/PatchedConicFixes/Plugins/`.

## License

MIT
