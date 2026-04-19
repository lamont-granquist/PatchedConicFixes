# PatchedConicFixes

This is a near overhaul of KSP 1.12.x's encounter resolution in the patched conic
solver.  It fixes several regressions introduced around KSP 1.8.0, along with
multiple design problems that KSP had from the start.  Many or most of the bugs
where SOI encounters to go undetected, even when your ship is on a clear intercept
course, have been fixed.  A few examples of affected scenarios include problems
with flickering orbits around Jool, and transfers to an ascending Moon in
RealSolarSystem, and many others. The replacement encounter solver is also
approximately 5× faster than the stock implementation.

## Known Limitations

This implementation has known failure conditions that stem from fundamental design
choices in the algorithm (shared with the old KSP implementation).  A meaningful fix
would require a ground-up rewrite (which is just dependent upon enough free time).

If you encounter still buggy behavior, it's quite likely a known issue. For the time being
bug reports for this mod are unnecessary.

## Requirements

- Kerbal Space Program 1.12.x
- KSPBurst (1.7.2 or compatible)
- Harmony (2.2.1 or compatible)

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
