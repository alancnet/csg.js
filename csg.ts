// # class Vector

// Represents a 3D vector.
// 
// Example usage:
// 
//     new CSG.Vector(1, 2, 3);
//     new CSG.Vector([1, 2, 3]);
//     new CSG.Vector({ x: 1, y: 2, z: 3 });
class Vector {

  x: number;
  y: number;
  z: number;

  constructor(x: number | number[] | Vector, y?: number, z?: number) {
    if (typeof x === 'number' && y && z) {
      this.x = x;
      this.y = y;
      this.z = z;
    } else if (x instanceof Vector) {
      this.x = x.x;
      this.y = x.y;
      this.z = x.z;
    } else if (x instanceof Array) {
      this.x = x[0];
      this.y = x[1];
      this.z = x[2];
    }
  }

  clone(): Vector {
    return new Vector(this.x, this.y, this.z);
  }

  negated(): Vector {
    return new Vector(-this.x, -this.y, -this.z);
  }

  plus(a: Vector): Vector {
    return new Vector(this.x + a.x, this.y + a.y, this.z + a.z);
  }

  minus(a: Vector): Vector {
    return new Vector(this.x - a.x, this.y - a.y, this.z - a.z);
  }

  times(a: number): Vector {
    return new Vector(this.x * a, this.y * a, this.z * a);
  }

  dividedBy(a: number): Vector {
    return new Vector(this.x / a, this.y / a, this.z / a);
  }

  dot(a: Vector): number {
    return this.x * a.x + this.y * a.y + this.z * a.z;
  }

  lerp(a: Vector, t: number): Vector {
    return this.plus(a.minus(this).times(t));
  }

  length(): number {
    return Math.sqrt(this.dot(this));
  }

  unit(): Vector {
    return this.dividedBy(this.length());
  }

  cross(a: Vector): Vector {
    return new Vector(
      this.y * a.z - this.z * a.y,
      this.z * a.x - this.x * a.z,
      this.x * a.y - this.y * a.x
    );
  }
}

// # class Vertex

// Represents a vertex of a polygon. Use your own vertex class instead of this
// one to provide additional features like texture coordinates and vertex
// colors. Custom vertex classes need to provide a `pos` property and `clone()`,
// `flip()`, and `interpolate()` methods that behave analogous to the ones
// defined by `CSG.Vertex`. This class provides `normal` so convenience
// functions like `CSG.sphere()` can return a smooth vertex normal, but `normal`
// is not used anywhere else.
class Vertex {
  constructor(public pos: Vector, public normal: Vector) {}

  clone(): Vertex {
    return new Vertex(this.pos.clone(), this.normal.clone());
  }

  // Invert all orientation-specific data (e.g. vertex normal). Called when the
  // orientation of a polygon is flipped.
  flip(): void {
    this.normal = this.normal.negated();
  }

  // Create a new vertex between this vertex and `other` by linearly
  // interpolating all properties using a parameter of `t`. Subclasses should
  // override this to interpolate additional properties.
  interpolate(other: Vertex, t: number): Vertex {
    return new Vertex(
      this.pos.lerp(other.pos, t),
      this.normal.lerp(other.normal, t)
    );
  }
};

// # class Plane

// Represents a plane in 3D space.
class Plane {
  constructor(public normal: Vector, public w: number) {}
  
  // `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
  // point is on the plane.
  static EPSILON = 1e-5;

  static fromPoints = function fromPoints(a: Vector, b: Vector, c: Vector): Plane {
    const n = b.minus(a).cross(c.minus(a)).unit();
    return new Plane(n, n.dot(a));
  }

  clone(): Plane {
    return new Plane(this.normal.clone(), this.w);
  }

  flip(): void {
    this.normal = this.normal.negated();
    this.w = -this.w;
  }

  // Split `polygon` by this plane if needed, then put the polygon or polygon
  // fragments in the appropriate lists. Coplanar polygons go into either
  // `coplanarFront` or `coplanarBack` depending on their orientation with
  // respect to this plane. Polygons in front or in back of this plane go into
  // either `front` or `back`.
  splitPolygon(polygon: Polygon, coplanarFront: Polygon[], coplanarBack: Polygon[], front: Polygon[], back: Polygon[]): void {
    const COPLANAR = 0;
    const FRONT = 1;
    const BACK = 2;
    const SPANNING = 3;

    // Classify each point as well as the entire polygon into one of the above
    // four classes.
    let polygonType = 0;
    const types = [];
    for (let i = 0; i < polygon.vertices.length; i++) {
      const t = this.normal.dot(polygon.vertices[i].pos) - this.w;
      const type = (t < -Plane.EPSILON) ? BACK : (t > Plane.EPSILON) ? FRONT : COPLANAR;
      polygonType |= type;
      types.push(type);
    }

    // Put the polygon in the correct list, splitting it when necessary.
    switch (polygonType) {
      case COPLANAR:
        (this.normal.dot(polygon.plane.normal) > 0 ? coplanarFront : coplanarBack).push(polygon);
        break;
      case FRONT:
        front.push(polygon);
        break;
      case BACK:
        back.push(polygon);
        break;
      case SPANNING:
        const f = [];
        const b = [];
        for (let i = 0; i < polygon.vertices.length; i++) {
          const j = (i + 1) % polygon.vertices.length;
          const ti = types[i];
          const tj = types[j];
          const vi = polygon.vertices[i];
          const vj = polygon.vertices[j];
          if (ti != BACK) f.push(vi);
          if (ti != FRONT) b.push(ti != BACK ? vi.clone() : vi);
          if ((ti | tj) == SPANNING) {
            const t = (this.w - this.normal.dot(vi.pos)) / this.normal.dot(vj.pos.minus(vi.pos));
            const v = vi.interpolate(vj, t);
            f.push(v);
            b.push(v.clone());
          }
        }
        if (f.length >= 3) {
          front.push(new Polygon(f, polygon.shared));
        }
        if (b.length >= 3) {
          back.push(new Polygon(b, polygon.shared));
        }
        break;
    }
  }
};

// # class Polygon

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
// instances but they must behave similarly (duck typing can be used for
// customization).
// 
// Each convex polygon has a `shared` property, which is shared between all
// polygons that are clones of each other or were split from the same polygon.
// This can be used to define per-polygon properties (such as surface color).
class Polygon {
  constructor(public vertices: Vertex[], public shared: any[] = []) {}

  plane = Plane.fromPoints(
    this.vertices[0].pos,
    this.vertices[1].pos,
    this.vertices[2].pos
  );

  clone() {
    const vertices = this.vertices.map(v => v.clone());
    return new Polygon(vertices, this.shared);
  }

  flip() {
    this.vertices.reverse().forEach(v => v.flip());
    this.plane.flip();
  }
}

// # class Node

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.
class CSGNode {

  plane: Plane;
  front: CSGNode;
  back: CSGNode;

  constructor(public polygons: Polygon[] = []) {
    this.build(polygons);
  }

  clone(): CSGNode {
    const node = new CSGNode();
    node.plane = this.plane && this.plane.clone();
    node.front = this.front && this.front.clone();
    node.back = this.back && this.back.clone();
    node.polygons = this.polygons.map(p => p.clone());
    return node;
  }

  // Convert solid space to empty space and empty space to solid space.
  invert(): void {
    for (let i = 0; i < this.polygons.length; i++) {
      this.polygons[i].flip();
    }
    this.plane.flip();
    if (this.front) this.front.invert();
    if (this.back) this.back.invert();
    const temp = this.front;
    this.front = this.back;
    this.back = temp;
  }

  // Recursively remove all polygons in `polygons` that are inside this BSP
  // tree.
  clipPolygons(polygons: Polygon[]): Polygon[] {
    if (!this.plane) {
      return polygons.slice();
    }
    let front = [];
    let back = [];
    for (let i = 0; i < polygons.length; i++) {
      this.plane.splitPolygon(polygons[i], front, back, front, back);
    }
    if (this.front) {
      front = this.front.clipPolygons(front);
    }
    if (this.back) {
      back = this.back.clipPolygons(back);
    } else {
      back = [];
    }
    return front.concat(back);
  }

  // Remove all polygons in this BSP tree that are inside the other BSP tree
  // `bsp`.
  clipTo(bsp: CSGNode): void {
    this.polygons = bsp.clipPolygons(this.polygons);
    if (this.front) {
      this.front.clipTo(bsp);
    }
    if (this.back) {
      this.back.clipTo(bsp);
    }
  }

  // Return a list of all polygons in this BSP tree.
  allPolygons(): Polygon[] {
    let polygons = this.polygons.slice();
    if (this.front) {
      polygons = polygons.concat(this.front.allPolygons());
    }
    if (this.back) {
      polygons = polygons.concat(this.back.allPolygons());
    }
    return polygons;
  }

  // Build a BSP tree out of `polygons`. When called on an existing tree, the
  // new polygons are filtered down to the bottom of the tree and become new
  // nodes there. Each set of polygons is partitioned using the first polygon
  // (no heuristic is used to pick a good split).
  build(polygons: Polygon[]): void {
    if (!polygons.length) {
      return;
    }

    if (!this.plane) {
      this.plane = polygons[0].plane.clone();
    }

    const front = [];
    const back = [];

    for (let i = 0; i < polygons.length; i++) {
      this.plane.splitPolygon(polygons[i], this.polygons, this.polygons, front, back);
    }

    if (front.length) {
      if (!this.front) this.front = new CSGNode();
      this.front.build(front);
    }

    if (back.length) {
      if (!this.back) this.back = new CSGNode();
      this.back.build(back);
    }
  }
}

interface CSGCylinderOptions {
  start?: Vector;
  end?: Vector;
  radius?: number;
  slices?: number;
}

interface CSGSphereOptions {
  center?: Vector;
  radius?: number;
  slices?: number;
  stacks?: number;
}

interface CSGCubeOptions {
  center?: Vector;
  radius?: number[] | number;
}

class CSG {
  constructor(private polygons: Polygon[] = []) {}

  // Construct an axis-aligned solid cuboid. Optional parameters are `center` and
  // `radius`, which default to `[0, 0, 0]` and `[1, 1, 1]`. The radius can be
  // specified using a single number or a list of three numbers, one for each axis.
  // 
  // Example code:
  // 
  //     var cube = CSG.cube({
  //       center: [0, 0, 0],
  //       radius: 1
  //     });
  static cube = function(options: CSGCubeOptions = {}): CSG {
    const c = new Vector(options.center || [0, 0, 0]);
    const r = !options.radius ? [1, 1, 1] : options.radius instanceof Array ?
             options.radius : [options.radius, options.radius, options.radius];
    
    return CSG.fromPolygons([
      [[0, 4, 6, 2], [-1, 0, 0]],
      [[1, 3, 7, 5], [+1, 0, 0]],
      [[0, 1, 5, 4], [0, -1, 0]],
      [[2, 6, 7, 3], [0, +1, 0]],
      [[0, 2, 3, 1], [0, 0, -1]],
      [[4, 5, 7, 6], [0, 0, +1]]
    ].map(info => {
      return new Polygon(info[0].map(function(i) {
        var pos = new Vector(
          c.x + r[0] * (2 * Number(!!(i & 1)) - 1),
          c.y + r[1] * (2 * Number(!!(i & 2)) - 1),
          c.z + r[2] * (2 * Number(!!(i & 4)) - 1)
        );
        return new Vertex(pos, new Vector(info[1]));
      }));
    }));
  }

  // Construct a solid sphere. Optional parameters are `center`, `radius`,
  // `slices`, and `stacks`, which default to `[0, 0, 0]`, `1`, `16`, and `8`.
  // The `slices` and `stacks` parameters control the tessellation along the
  // longitude and latitude directions.
  // 
  // Example usage:
  // 
  //     var sphere = CSG.sphere({
  //       center: [0, 0, 0],
  //       radius: 1,
  //       slices: 16,
  //       stacks: 8
  //     });
  static sphere = function(options: CSGSphereOptions = {}): CSG {
    const c = new Vector(options.center || [0, 0, 0]);
    const r = options.radius || 1;
    const slices = options.slices || 16;
    const stacks = options.stacks || 8;
    const polygons: Polygon[] = [];
    let vertices: Vertex[] = [];
    
    function vertex(theta, phi) {
      const dir = new Vector(
        Math.cos(theta * Math.PI * 2) * Math.sin(phi * Math.PI),
        Math.cos(phi * Math.PI),
        Math.sin(theta * Math.PI * 2) * Math.sin(phi * Math.PI)
      );
      vertices.push(new Vertex(c.plus(dir.times(r)), dir));
    }
  
    for (let i = 0; i < slices; i++) {
      for (let j = 0; j < stacks; j++) {
        vertices = [];
        vertex(i / slices, j / stacks);
        if (j > 0) vertex((i + 1) / slices, j / stacks);
        if (j < stacks - 1) vertex((i + 1) / slices, (j + 1) / stacks);
        vertex(i / slices, (j + 1) / stacks);
        polygons.push(new Polygon(vertices));
      }
    }

    return CSG.fromPolygons(polygons);
  }

  // Construct a solid cylinder. Optional parameters are `start`, `end`,
  // `radius`, and `slices`, which default to `[0, -1, 0]`, `[0, 1, 0]`, `1`, and
  // `16`. The `slices` parameter controls the tessellation.
  // 
  // Example usage:
  // 
  //     var cylinder = CSG.cylinder({
  //       start: [0, -1, 0],
  //       end: [0, 1, 0],
  //       radius: 1,
  //       slices: 16
  //     });
  static cylinder = function(options: CSGCylinderOptions = {}): CSG {
    const s = new Vector(options.start || [0, -1, 0]);
    const e = new Vector(options.end || [0, 1, 0]);
    const ray = e.minus(s);
    const r = options.radius || 1;
    const slices = options.slices || 16;
    const axisZ = ray.unit(), isY = (Math.abs(axisZ.y) > 0.5);
    const axisX = new Vector(Number(isY), Number(!isY), 0).cross(axisZ).unit();
    const axisY = axisX.cross(axisZ).unit();
    const start = new Vertex(s, axisZ.negated());
    const end = new Vertex(e, axisZ.unit());
    const polygons: Polygon[] = [];

    function point(stack: number, slice: number, normalBlend: number): Vertex {
      const angle = slice * Math.PI * 2;
      const out = axisX.times(Math.cos(angle)).plus(axisY.times(Math.sin(angle)));
      const pos = s.plus(ray.times(stack)).plus(out.times(r));
      const normal = out.times(1 - Math.abs(normalBlend)).plus(axisZ.times(normalBlend));
      return new Vertex(pos, normal);
    }

    for (let i = 0; i < slices; i++) {
      const t0 = i / slices, t1 = (i + 1) / slices;
      polygons.push(new Polygon([start, point(0, t0, -1), point(0, t1, -1)]));
      polygons.push(new Polygon([point(0, t1, 0), point(0, t0, 0), point(1, t0, 0), point(1, t1, 0)]));
      polygons.push(new Polygon([end, point(1, t1, 1), point(1, t0, 1)]));
    }

    return CSG.fromPolygons(polygons);
  }

  // Construct a CSG solid from a list of `CSG.Polygon` instances.
  static fromPolygons = function(polygons: Polygon[]): CSG {
    return new CSG(polygons);
  }

  static Node = CSGNode;
  static Polygon = Polygon;
  static Plane = Plane;
  static Vertex = Vertex;
  static Vector = Vector;

  clone(): CSG {
    return new CSG(this.polygons.map(p => p.clone()));
  }

  toPolygons(): Polygon[] {
    return this.polygons;
  }

  // Return a new CSG solid representing space in either this solid or in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.union(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |       +----+
  //     +----+--+    |       +----+       |
  //          |   B   |            |       |
  //          |       |            |       |
  //          +-------+            +-------+
  // 
  union(csg: CSG): CSG {
    const a = new CSGNode(this.clone().polygons);
    const b = new CSGNode(csg.clone().polygons);
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.build(b.allPolygons());
    return CSG.fromPolygons(a.allPolygons());
  }

  // Return a new CSG solid representing space in this solid but not in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.subtract(B)
  // 
  //     +-------+            +-------+
  //     |       |            |       |
  //     |   A   |            |       |
  //     |    +--+----+   =   |    +--+
  //     +----+--+    |       +----+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  subtract(csg: CSG): CSG {
    const a = new CSGNode(this.clone().polygons);
    const b = new CSGNode(csg.clone().polygons);
    a.invert();
    a.clipTo(b);
    b.clipTo(a);
    b.invert();
    b.clipTo(a);
    b.invert();
    a.build(b.allPolygons());
    a.invert();
    return CSG.fromPolygons(a.allPolygons());
  }

  // Return a new CSG solid representing space both this solid and in the
  // solid `csg`. Neither this solid nor the solid `csg` are modified.
  // 
  //     A.intersect(B)
  // 
  //     +-------+
  //     |       |
  //     |   A   |
  //     |    +--+----+   =   +--+
  //     +----+--+    |       +--+
  //          |   B   |
  //          |       |
  //          +-------+
  // 
  intersect(csg: CSG): CSG {
    var a = new CSGNode(this.clone().polygons);
    var b = new CSGNode(csg.clone().polygons);
    a.invert();
    b.clipTo(a);
    b.invert();
    a.clipTo(b);
    b.clipTo(a);
    a.build(b.allPolygons());
    a.invert();
    return CSG.fromPolygons(a.allPolygons());
  }

  // Return a new CSG solid with solid and empty space switched. This solid is
  // not modified.
  inverse(): CSG {
    const csg = this.clone();
    csg.polygons.forEach(p => p.flip());
    return csg;
  }
}
